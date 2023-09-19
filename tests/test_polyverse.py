import pytest
import pandas as pd
from copy import deepcopy

from polyverse.utils import canonical_smiles
from polyverse import constants as ks
from polyverse import polymerize

npartitions = 1


@pytest.fixture
def example_data():
    smiles_dict = {
        "adipic acid": "C(CCC(=O)O)CC(=O)O",
        "hexane diamine": "C(CCCN)CCN",
        "3-Aminobenzoyl chloride": "C1=CC(=CC(=C1)N)C(=O)Cl",
        "Bicyclo[2.2.1]heptane-2,5-dicarboxylic acid": "O=C(O)C1CC2CC1CC2C(=O)O",
        "Cyclohexene": "C1=CCCCC1",
        "Cyclopentadiene": "C1=CCC=C1",
        "TFTPN": canonical_smiles(
            "C(#N)C1=C(C(=C(C(=C1F)F)C#N)F)F"
        ),  # Tetrafluoroterephthalonitrile
        "TCTPN": canonical_smiles(
            "C(#N)C1=C(C(=C(C(=C1Cl)Cl)C#N)Cl)Cl"  # Tetrachloroterephthalonitrile
        ),
        "TTSBI": canonical_smiles(
            "CC1(CC2(CC(C3=CC(=C(C=C32)O)O)(C)C)C4=CC(=C(C=C41)O)O)C"
        ),  # 5,5',6,6'-TETRAHYDROXY-3,3,3',3'-TETRAMETHYL-1,1'-SPIROBISINDANE
        "3,4,5-Trihydroxybenzhydrazide": "NNC(=O)C1=CC(=C(C(=C1)O)O)O",
        "3,4,5-Trifluorobenzaldehyde": "O=Cc1cc(F)c(F)c(F)c1",
        "maleimide": "C1=CC(=O)NC1=O",
        "3-Chlorocyclopentene": "C1CC(C=C1)Cl",
        "2,3,5,6-tetrafluoro bromobenzene": "C1=C(C(=C(Br)C(=C1F)F)F)F",
        "allylethylene": canonical_smiles("C=CCC=C"),
    }
    # DO NOT DELETE ITEMS FROM id_sm. ONLY ADD.
    id_sm = [
        ("ZINC0", smiles_dict["adipic acid"]),
        ("ZINC1", smiles_dict["hexane diamine"]),
        ("ZINC2", "C(C(N=C=O)CCN)CCN"),
        ("ZINC3", smiles_dict["3-Aminobenzoyl chloride"]),
        ("ZINC4", "C12C(CC(CC1CN)C2)C(=O)O"),  # monomer with norbornane
        ("ZINC5", smiles_dict["Bicyclo[2.2.1]heptane-2,5-dicarboxylic acid"]),
        (
            "ZINC6",
            canonical_smiles("O=C1C2C3C=CC(O3)C2C(=O)N1c1cc(C(F)(F)F)cc(C(F)(F)F)c1"),
        ),
        ("ZINC7", smiles_dict["Cyclohexene"]),
        ("ZINC8", smiles_dict["Cyclopentadiene"]),
        ("ZINC9", smiles_dict["TFTPN"]),
        ("ZINC10", smiles_dict["TTSBI"]),
        ("ZINC11", smiles_dict["3,4,5-Trihydroxybenzhydrazide"]),
        ("ZINC12", smiles_dict["3,4,5-Trifluorobenzaldehyde"]),
        ("ZINC13", smiles_dict["maleimide"]),
        ("ZINC14", smiles_dict["3-Chlorocyclopentene"]),
        ("ZINC15", smiles_dict["2,3,5,6-tetrafluoro bromobenzene"]),
        ("ZINC16", "C1=C(C(=C([OH])C(=C1F)F)F)F"),
        ("ZINC17", smiles_dict["TCTPN"]),
        (
            "ZINC18",
            canonical_smiles(
                "C(#N)C1=C(C(=C(C(=C1Br)Br)C#N)Br)Br"  # same as TCTPN, except sub Cl with Br.
            ),
        ),
        (
            "ZINC19",
            canonical_smiles(
                "C(#N)C1=C(C(=C(C(=C1I)I)C#N)I)I"  # same as TCTPN, except sub Cl with I.
            ),
        ),
        (
            "ZINC20",
            canonical_smiles(
                "O=C(O)CNC(=O)C5=C(O)C3=C(Cl)C(=C(Cl)[N]3N(CC4=CC=CC=C4)C5=O)Cl"
            ),
        ),
        (
            "ZINC21",
            canonical_smiles("Clc1sc(I)c(Cl)c1Cl"),
        ),
        (
            "ZINC22",
            canonical_smiles("FC1=C(F)C(Cl)=C(F)C(F)=C1Cl"),
        ),
        (
            "ZINC23",
            canonical_smiles(
                "N#Cc3c(F)c(F)c(C#N)c(Nc2ccc(Nc1c(F)c(C#N)c(F)c(F)c1C#N)cc2)c3F"
            ),
        ),
        (
            "ZINC24",
            canonical_smiles("CC4(C)CC2(C[C+](C)(C)(C)c1cc(O)c(O)cc12)c3cc(O)c(O)cc34"),
        ),
        (
            "ZINC25",
            smiles_dict["allylethylene"],
        ),
        ("ZINC26", canonical_smiles("C=CCCC(C)CCC=C")),
        (
            "ZINC27",
            canonical_smiles("C=CCCC(C)CC(O)C=C"),
        ),
        ("ZINC28", canonical_smiles("C#CCC#C")),
        ("ZINC29", canonical_smiles("[N-]=[N+]=NCOCOCN=[N+]=[N-]")),
    ]
    index = list(range(len(id_sm)))
    df = pd.DataFrame(
        {
            "index": index,
            "id-smiles": id_sm,
        }
    )

    return {"df": df, "smiles_dict": smiles_dict}


def test_romp(example_data):
    correct_df = pd.DataFrame(
        {
            "zinc_ids": ["ZINC6"],
            "polymer": ["*=CC1OC(C=*)C2C(=O)N(c3cc(C(F)(F)F)cc(C(F)(F)F)c3)C(=O)C12"],
            "reactants": [
                canonical_smiles(
                    "O=C3C2C1C=CC(O1)C2C(=O)N3c4cc(C(F)(F)F)cc(C(F)(F)F)c4"
                ),
            ],
        }
    )
    molecules = example_data["df"].iloc[[6, 7, 8, 13, 14], :]
    romp = polymerize.ROMP(npartitions=npartitions)
    polymers = romp.generate_polymers(molecules)
    # check that we get the polymers we expect
    pd.testing.assert_frame_equal(polymers, correct_df)
