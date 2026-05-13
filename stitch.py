import json
import re
import argparse
from rdkit import Chem

# =========================================================
# NATURAL AA DATABASE
# =========================================================
NATURAL = {
    "A": "C[C@H](N)C(O)=O",
    "G": "NCC(O)=O",
    "C": "N[C@@H](CS)C(O)=O",
    "D": "N[C@@H](CC(O)=O)C(O)=O",
    "E": "N[C@@H](CCC(O)=O)C(O)=O",
    "F": "N[C@@H](Cc1ccccc1)C(O)=O",
    "H": "N[C@@H](Cc1cnc[nH]1)C(O)=O",
    "I": "CC[C@H](C)[C@H](N)C(O)=O",
    "K": "N[C@@H](CCCCN)C(O)=O",
    "L": "CC(C)C[C@H](N)C(O)=O",
    "M": "CSCC[C@H](N)C(O)=O",
    "N": "N[C@@H](CC(N)=O)C(O)=O",
    "P": "OC(=O)[C@@H]1CCCN1",
    "Q": "N[C@@H](CCC(N)=O)C(O)=O",
    "R": "N[C@@H](CCCNC(N)=N)C(O)=O",
    "S": "N[C@@H](CO)C(O)=O",
    "T": "C[C@@H](O)[C@H](N)C(O)=O",
    "V": "CC(C)[C@H](N)C(O)=O",
    "W": "N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O",
    "Y": "N[C@@H](Cc1ccc(O)cc1)C(O)=O",
}

# =========================================================
# BACKBONE DETECTION
# =========================================================

def find_backbone_n_idx(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        for alpha in atom.GetNeighbors():
            if alpha.GetSymbol() != "C":
                continue
            for carbonyl in alpha.GetNeighbors():
                if carbonyl.GetSymbol() != "C":
                    continue
                for b in carbonyl.GetBonds():
                    other = b.GetOtherAtom(carbonyl)
                    if other.GetSymbol() == "O" and b.GetBondType() == Chem.BondType.DOUBLE:
                        return atom.GetIdx()

    nitrogens = [a for a in mol.GetAtoms() if a.GetSymbol() == "N"]
    carboxyl_cs = _carboxyl_carbon_indices(mol)
    if not nitrogens or not carboxyl_cs:
        raise RuntimeError("find_backbone_n: no backbone N found")

    best_idx, min_dist = None, 999
    for n in nitrogens:
        for c_idx in carboxyl_cs:
            path = Chem.GetShortestPath(mol, n.GetIdx(), c_idx)
            if path and len(path) - 1 < min_dist:
                min_dist = len(path) - 1
                best_idx = n.GetIdx()
    if best_idx is None:
        raise RuntimeError("find_backbone_n: no backbone N found")
    return best_idx


def _carboxyl_carbon_indices(mol):
    result = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        has_double_o = has_single_o = False
        for b in atom.GetBonds():
            o = b.GetOtherAtom(atom)
            if o.GetSymbol() != "O":
                continue
            if b.GetBondType() == Chem.BondType.DOUBLE:
                has_double_o = True
            elif b.GetBondType() == Chem.BondType.SINGLE:
                has_single_o = True
        if has_double_o and has_single_o:
            result.append(atom.GetIdx())
    return result


def find_backbone_c_idx(mol):
    candidates = _carboxyl_carbon_indices(mol)
    if not candidates:
        raise RuntimeError("find_backbone_c: no carboxyl C found")
    if len(candidates) == 1:
        return candidates[0]
    n_idx = find_backbone_n_idx(mol)
    best_idx, min_dist = None, 999
    for c_idx in candidates:
        path = Chem.GetShortestPath(mol, n_idx, c_idx)
        if path and len(path) - 1 < min_dist:
            min_dist = len(path) - 1
            best_idx = c_idx
    return best_idx


# =========================================================
# TAGGING & STITCHING
# =========================================================

def tag_n_terminus(mol):
    rw = Chem.RWMol(mol)
    n_idx = find_backbone_n_idx(rw)
    star = Chem.Atom("*")
    star.SetAtomMapNum(1)
    star_idx = rw.AddAtom(star)
    rw.AddBond(n_idx, star_idx, Chem.BondType.SINGLE)
    return rw


def tag_c_terminus(mol):
    rw = Chem.RWMol(mol)
    c_idx = find_backbone_c_idx(rw)
    oh_idx = None
    for b in rw.GetAtomWithIdx(c_idx).GetBonds():
        o = b.GetOtherAtom(rw.GetAtomWithIdx(c_idx))
        if o.GetSymbol() == "O" and b.GetBondType() == Chem.BondType.SINGLE:
            oh_idx = o.GetIdx()
            break
    if oh_idx is not None:
        rw.RemoveAtom(oh_idx)
        if oh_idx < c_idx:
            c_idx -= 1
    star = Chem.Atom("*")
    star.SetAtomMapNum(2)
    star_idx = rw.AddAtom(star)
    rw.AddBond(c_idx, star_idx, Chem.BondType.SINGLE)
    return rw


def prepare_residue(smiles, res_idx, tag_n=True, tag_c=True):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    for atom in mol.GetAtoms():
        atom.SetIntProp("res_idx", res_idx)
    if tag_n:
        mol = tag_n_terminus(mol)
    if tag_c:
        mol = tag_c_terminus(mol)
    return mol


def stitch(left, right):
    combo = Chem.CombineMols(left, right)
    rw = Chem.RWMol(combo)
    star_atoms = [(a.GetAtomMapNum(), a.GetIdx(), a.GetNeighbors()[0].GetIdx())
                  for a in rw.GetAtoms() if a.GetSymbol() == "*"]
    c_map  = next(s for s in star_atoms if s[0] == 2)
    n_map  = next(s for s in star_atoms if s[0] == 1)
    c_star_idx, c_atom_idx = c_map[1], c_map[2]
    n_star_idx, n_atom_idx = n_map[1], n_map[2]
    rw.AddBond(c_atom_idx, n_atom_idx, Chem.BondType.SINGLE)
    hi, lo = (c_star_idx, n_star_idx) if c_star_idx > n_star_idx else (n_star_idx, c_star_idx)
    rw.RemoveAtom(hi)
    rw.RemoveAtom(lo)
    return rw.GetMol()


# =========================================================
# PARSERS & CYCLIZATION
# =========================================================

def parse_mod(seq, i):
    j = seq.index("}", i)
    content = seq[i + 1:j]
    if ":" in content:
        kind, code = content.split(":", 1)
    else:
        kind = code = content
    return kind.lower(), code.lower(), j + 1


def parse_sequence(seq):
    residues, nt, ct, i = [], None, None, 0
    while i < len(seq):
        if seq[i] == "{":
            kind, code, i = parse_mod(seq, i)
            if kind == "nt": nt = code
            elif kind == "ct": ct = code
            elif kind == "ptm":
                aa, _, is_d = residues[-1]
                residues[-1] = (aa, code, is_d)
            elif kind == "nnr": residues.append(("NNR", code.upper(), False))
            elif kind == "d":
                aa, mod, _ = residues[-1]
                residues[-1] = (aa, mod, True)
            elif kind == "cyc": pass
            else: raise ValueError(f"Unknown modification type: '{kind}'")
            continue
        residues.append((seq[i].upper(), None, False))
        i += 1
    return residues, nt, ct


def parse_cyclization(seq):
    match = re.search(r"\{cyc:([^}]*)\}", seq.lower())
    if not match: return None
    content = match.group(1)
    if content == "n-c": return {"type": "head_tail"}
    pairs = []
    for part in content.split(","):
        a, b = part.split("-")
        pairs.append((int(a) - 1, int(b) - 1))
    return {"type": "pairs", "pairs": pairs}


def apply_cyclization(mol, cyc_info):
    if not cyc_info:
        return mol

    rw = Chem.RWMol(mol)

    if cyc_info["type"] == "head_tail":

        atoms = list(rw.GetAtoms())
        n_atom_idx = None

        # Find N-terminus backbone nitrogen
        for a in atoms:
            if a.GetSymbol() != "N" or a.GetIntProp("res_idx") != 0:
                continue

            for alpha in a.GetNeighbors():
                if alpha.GetSymbol() != "C":
                    continue

                for carbonyl in alpha.GetNeighbors():
                    if carbonyl.GetSymbol() != "C":
                        continue

                    for b in carbonyl.GetBonds():
                        other = b.GetOtherAtom(carbonyl)

                        if (
                            other.GetSymbol() == "O"
                            and b.GetBondType() == Chem.BondType.DOUBLE
                        ):
                            n_atom_idx = a.GetIdx()
                            break

                    if n_atom_idx is not None:
                        break

                if n_atom_idx is not None:
                    break

        # Fallback for unusual residues
        if n_atom_idx is None:
            for a in atoms:
                if (
                    a.GetSymbol() == "N"
                    and a.GetIntProp("res_idx") == 0
                    and not a.IsInRing()
                ):
                    n_atom_idx = a.GetIdx()
                    break

        if n_atom_idx is None:
            raise RuntimeError(
                "Head-tail cyclization: could not find backbone N."
            )

        # Find C-terminus carbonyl carbon
        max_res = max(
            a.GetIntProp("res_idx")
            for a in atoms
            if a.HasProp("res_idx")
        )

        c_atom_idx = None

        for a in atoms:
            if (
                a.GetSymbol() == "C"
                and a.GetIntProp("res_idx") == max_res
            ):
                for b in a.GetBonds():
                    o = b.GetOtherAtom(a)

                    if (
                        o.GetSymbol() == "O"
                        and b.GetBondType() == Chem.BondType.DOUBLE
                    ):
                        c_atom_idx = a.GetIdx()

        if c_atom_idx is None:
            raise RuntimeError(
                "Head-tail cyclization: could not find C-terminal C=O."
            )

        # Remove terminal OH
        oh_idx = None

        for b in rw.GetAtomWithIdx(c_atom_idx).GetBonds():
            o = b.GetOtherAtom(rw.GetAtomWithIdx(c_atom_idx))

            if (
                o.GetSymbol() == "O"
                and b.GetBondType() == Chem.BondType.SINGLE
            ):
                oh_idx = o.GetIdx()
                break

        if oh_idx is not None:
            rw.RemoveAtom(oh_idx)

            if oh_idx < c_atom_idx:
                c_atom_idx -= 1

            if oh_idx < n_atom_idx:
                n_atom_idx -= 1

        rw.AddBond(
            n_atom_idx,
            c_atom_idx,
            Chem.BondType.SINGLE
        )

    elif cyc_info["type"] == "pairs":

        sulfurs = {}

        # Collect ONLY cysteine-like thiol sulfurs
        for a in rw.GetAtoms():

            if a.GetSymbol() != "S":
                continue

            if not a.HasProp("res_idx"):
                continue

            # Cysteine sulfur has ONE carbon neighbor
            # Methionine sulfur has TWO carbon neighbors
            carbon_neighbors = [
                n for n in a.GetNeighbors()
                if n.GetSymbol() == "C"
            ]

            if len(carbon_neighbors) != 1:
                continue

            sulfurs.setdefault(
                a.GetIntProp("res_idx"),
                []
            ).append(a.GetIdx())

        # Create disulfide bonds
        for r1, r2 in cyc_info["pairs"]:

            # Validate residue 1
            if r1 not in sulfurs:
                raise ValueError(
                    f"Cyclization error: residue {r1 + 1} "
                    f"is not a cysteine."
                )

            # Validate residue 2
            if r2 not in sulfurs:
                raise ValueError(
                    f"Cyclization error: residue {r2 + 1} "
                    f"is not a cysteine."
                )

            rw.AddBond(
                sulfurs[r1][0],
                sulfurs[r2][0],
                Chem.BondType.SINGLE
            )

    return rw.GetMol()


# =========================================================
# DATABASE
# =========================================================

def extract_aa_letter(field):
    return field.split("/")[-1].strip()

def load_databases(path):
    with open(path) as f: data = json.load(f)
    user_db, map_db, nnr_db = {}, {}, {}
    for e in data:
        smiles = e["SMILES"]
        aa = extract_aa_letter(e["Natural Amino Acid"]) if e.get("Natural Amino Acid") else None
        for key in ("MAP notation", "User Input"):
            tag = e.get(key, "")
            if tag and aa and "{" in tag:
                inside = tag[tag.index("{") + 1: tag.index("}")]
                if ":" in inside:
                    kind, code = inside.split(":", 1)
                    map_db[(aa, kind.lower(), code.lower())] = smiles
        if e.get("User Code") and aa:
            user_db[(aa, e["User Code"].lower())] = smiles
        if e.get("MAP notation", "").startswith("{nnr:"):
            content = e["MAP notation"].strip("{}")
            if ":" in content:
                _, code = content.split(":", 1)
                nnr_db[code.upper()] = smiles
    return user_db, map_db, nnr_db

def resolve_mod(aa, kind, code, user_db, map_db):
    key = (aa, kind.lower(), code.lower())
    if key in map_db: return map_db[key]
    if kind.lower() == "ptm" and (aa, code.lower()) in user_db:
        return user_db[(aa, code.lower())]
    raise KeyError(f"No modification found for {aa}:{kind}:{code}")


# =========================================================
# MAIN ENGINE
# =========================================================

def sequence_to_smiles(seq, user_db, map_db, nnr_db):
    residues, nt, ct = parse_sequence(seq)
    cyc_info = parse_cyclization(seq)
    mols = []

    for i, (aa, mod, is_d) in enumerate(residues):
        if aa == "NNR": smiles = nnr_db[mod]
        elif mod: smiles = resolve_mod(aa, "ptm", mod, user_db, map_db)
        else: smiles = NATURAL[aa]

        if i == 0 and nt: smiles = resolve_mod(aa, "nt", nt, user_db, map_db)
        elif i == len(residues) - 1 and ct: smiles = resolve_mod(aa, "ct", ct, user_db, map_db)

        mols.append(prepare_residue(smiles, i, tag_n=(i != 0), tag_c=(i != len(residues) - 1)))

    peptide = mols[0]
    for nxt in mols[1:]:
        peptide = stitch(peptide, nxt)

    peptide = apply_cyclization(peptide, cyc_info)

    # FIX: call UpdatePropertyCache at the molecule level, not in a per-atom loop.
    # A per-atom loop processes each nitrogen without the context of the full
    # aromatic ring system, causing RDKit to strip [nH] from histidine's imidazole
    # (and any other aromatic NH). The molecule-level call preserves it correctly.
    peptide.UpdatePropertyCache(strict=False)

    Chem.FastFindRings(peptide)
    Chem.SanitizeMol(
        peptide,
        Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
    )

    Chem.AssignStereochemistry(peptide, cleanIt=True, force=True)
    for atom in peptide.GetAtoms():
        if not atom.HasProp("res_idx"): continue
        rid = atom.GetIntProp("res_idx")
        if atom.GetSymbol() != "C" or atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED: continue
        if rid < len(residues) and residues[rid][2]:
            flipped = (Chem.ChiralType.CHI_TETRAHEDRAL_CCW if atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW
                       else Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            atom.SetChiralTag(flipped)

    Chem.AssignStereochemistry(peptide, cleanIt=True, force=True)
    return Chem.MolToSmiles(peptide, isomericSmiles=True)


def read_sequences(file_path):
    sequences = []

    if file_path.lower().endswith((".fasta", ".fa")):
        with open(file_path, encoding="utf-8-sig") as f:
            header = None
            seq = ""

            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    # save previous entry
                    if seq:
                        sequences.append((header, seq))

                    header = line[1:]   # remove ">"
                    seq = ""

                else:
                    seq += line

            # save last entry
            if seq:
                sequences.append((header, seq))

    else:
        with open(file_path, encoding="utf-8-sig") as f:
            for i, line in enumerate(f, start=1):
                s = line.strip()
                if s:
                    # no header in txt files
                    sequences.append((f"seq_{i}", s))

    return sequences


# =========================================================
# CLI
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="Convert peptide sequences to isomeric SMILES."
    )
    parser.add_argument("sequence", nargs="?", help="Single sequence string")
    parser.add_argument("--db", default="data/merged_final_data.json",
                        help="Path to modifications JSON database")
    parser.add_argument("-i", "--input",  help="Input file (.txt or .fasta)")
    parser.add_argument("-o", "--output", help="Output CSV file")
    args = parser.parse_args()

    user_db, map_db, nnr_db = load_databases(args.db)

    if args.input and args.output:
        sequences = read_sequences(args.input)

        with open(args.output, "w") as fout:
            fout.write("ID,sequence,SMILES\n")

            errors = []

            for header, seq in sequences:
                try:
                    smi = sequence_to_smiles(seq, user_db, map_db, nnr_db)
                    fout.write(f"{header},{seq},{smi}\n")

                except Exception as e:
                    errors.append((header, seq, str(e)))

            if errors:
                err_file = args.output + ".errors"

                with open(err_file, "w") as ef:
                    ef.write("header,sequence,error\n")

                    for header, seq, err in errors:
                        ef.write(f"{header},{seq},{err}\n")

                print(f"Skipped {len(errors)} sequences — see {err_file}")

        print(f"Done. Written to {args.output}")

    elif args.sequence:
        smi = sequence_to_smiles(args.sequence, user_db, map_db, nnr_db)
        print(smi)

    else:
        parser.error("Provide a sequence or use -i / -o for batch mode.")


if __name__ == "__main__":
    main()