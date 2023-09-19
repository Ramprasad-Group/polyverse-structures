from rdkit import Chem


def canonicalSmiles(sm):
    return Chem.MolToSmiles(Chem.MolFromSmiles(sm))
