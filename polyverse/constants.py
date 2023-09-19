romp_rxn_smarts_sub = "&!$(*C(=O)*)&!$(*C[F,Cl,Br,I])&!$(*C(=O)[OH])&!$(*CO[F,Cl,Br,I])"
romp_rxn_smarts = (
    f"[CH1R{romp_rxn_smarts_sub}:0]="
    + f"[CH1R{romp_rxn_smarts_sub}:1]>>([#0]=[C:0].[#0]=[C:1])"
)
