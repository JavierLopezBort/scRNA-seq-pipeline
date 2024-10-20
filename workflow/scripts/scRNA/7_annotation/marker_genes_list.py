Persio = {
    "SPG": ["MAGEA4","DMRT1","SOX4"],
    "SCT_lep": ["DPH7","SYCP3","SCML1","TEX19"],
    "SCT_zyg": ["LY6K","SELENOT","TDRG1"],
    "SCT_pach": ["PIWIL1","CCDC112"],
    "SCT_dyp": ["AURKA","CCNA1"],
    "STD_early": ["TEX29"],
    "STD_late": ["PRM2","TNP1","PRM1"],
    "PMCs": ["ACTA2","MYH11"],
    "Fibrotic_PMCs": ["CFD","CLEC3B"],
    "Sertoli_cells": ["FATE1","CITED1"],
    "Leydig_cells": ["HSD17B3","STAR","INSL3"],
    "Endothelial_cells": ["VWF"],
    "Macrophages": ["LYZ","C1QA","CD14"]
}

Persio_meiosis = {
    "SPG": "meiotic",
    "SCT_lep": "meiotic",
    "SCT_zyg": "meiotic",
    "SCT_pach": "meiotic",
    "SCT_dyp": "meiotic",
    "STD_early": "meiotic",
    "STD_late": "meiotic",
    "PMCs": "somatic",
    "Fibrotic_PMCs": "somatic",
    "Sertoli_cells": "somatic",
    "Leydig_cells": "somatic",
    "Endothelial_cells": "somatic",
    "Macrophages": "somatic"
}

Hermann = {
    "SSC": ["ID4","NANOS2","TCN2"],
    "SPG_prog": ["EGR4","PLPPR5","TSPAN33"],
    "SPG_early_diff": ["ASB9","NANOS3","TMEM55A"],
    "SPG_late_diff": ["CDT1","ISOC1","NMT2"],
    "SCT_prelep": ["CT55","DHRS13","HSD17B14"],
    "SCT_lep_zyg": ["HIST3H3","MEIOB","TEX101"],
    "SCT_pach": ["C9orf57","CETN3","MGAT4D"],
    "SCT_dyp": ["ARMC1","PHOSPHO2","TPPP3"],
    "STD_early_r": ["C17orf98","ENPP2","LRRC3B"],
    "STD_mid_r": ["PRSS37","PRSS58","TP53TG5"],
    "STD_late_r": ["C17orf74","FSCN3","PHOSPHO1"],
    "PMCs": ["ACTA2","MYH11"],
    "Fibrotic_PMCs": ["CFD","CLEC3B"],
    "Sertoli_cells": ["FATE1","CITED1"],
    "Leydig_cells": ["HSD17B3","STAR","INSL3"],
    "Endothelial_cells": ["VWF"],
    "Macrophages": ["LYZ","C1QA","CD14"]
}

Hermann_meiosis = {
    "SSC": "meiotic",
    "SPG_prog": "meiotic",
    "SPG_early_diff": "meiotic",
    "SPG_late_diff": "meiotic",
    "SCT_prelep": "meiotic",
    "SCT_lep_zyg": "meiotic",
    "SCT_pach": "meiotic",
    "SCT_dyp": "meiotic",
    "STD_early_r": "meiotic",
    "STD_mid_r": "meiotic",
    "STD_late_r": "meiotic",
    "PMCs": "somatic",
    "Fibrotic_PMCs": "somatic",
    "Sertoli_cells": "somatic",
    "Leydig_cells": "somatic",
    "Endothelial_cells": "somatic",
    "Macrophages": "somatic"
}

Murat = {
    "Sertoli_cells": ["CLU"],
    "PSMCs": ["TAGLN", "ACTA2"],
    "Endothelial cells": ["CD34", "TM4SF1"],
    "Macrophages": ["APOE", "CD74"],
    "Leydig_cells": ["STAR", "CYP11A1"],
    "SPG_undiff": ["GFRA1", "PIWIL4", "STRA8"],
    "SPG_diff": ["DMRT1", "STRA8"],
    "SCT_lep": ["SYCE1", "SYCP2", "TANK", "AURKA"],
    "SCT_zyg": ["SYCP1", "SYCP2", "TANK", "AURKA"],
    "SCT_pach": ["PIWIL1", "SYCP2", "TANK", "AURKA"],
    "STD_early_r": ["LRRIQ1"],
    "STD_late_r": ["ACRV1", "SPACA1"],
    "STD_e": ["SPATA3", "NRBP1", "PRM1", "GABBR2"]
}

Murat_meiosis = {
    "Sertoli_cells": "somatic",
    "PSMCs": "somatic",
    "Endothelial cells": "somatic",
    "Macrophages": "somatic",
    "Leydig_cells": "somatic",
    "SPG_undiff": "meiotic",
    "SPG_diff": "meiotic",
    "SCT_lep": "meiotic",
    "SCT_zyg": "meiotic",
    "SCT_pach": "meiotic",
    "STD_early_r": "meiotic",
    "STD_late_r": "meiotic",
    "STD_e": "meiotic"
}

Taelman = {
    "Epithelial rete testis": ["PAX8","KRT19", "PCP4"],
    "Sertoli": ["AMH", "SOX9", "GATM", "DMRT1"],
    "Endothelial cells": ["PECAM1"],
    "Smooth muscle cells": ["RGS5"],
	"Immune cells": ["CD53"],
    "Epithelial cells from mesonephroi/epididymis": ["EPCAM", "ANXA2"],
	"Stromal cells from mesonephroi/epididymis": ["NR2F2", "DLK1", "PDGFRA"],
	"Stromal cells from the testis": ["NR2F2","DLK1","PDGFRA","ARX","GATA4"],
	"Leydig cells": ["CYP17A1", "INSL3"],
	"Epithelial cells of the MTs and WD": ["LHX1", "PAX2"],
	"Secretory tubules": ["AQP3", "CALB2", "GJA1", "CDH2", "AMH"],
	"Proliferating Sertoli cells, gonadal stromal cells, and mesonephric/epididymal stromal cells": ["MKI67"],
	"Coelomic/ovarian surface epithelium": ["UPK3B", "KRT19", "WT1"],
	"Pre-granulosa cells": ["FOXL2", "GATM", "WT1"],
	"Tubular cells of the f.MTs, f.MD, and f.WD": ["EPCAM", "LHX1", "PAX2"],
	"Mesonephric podocytes": ["PODXL","CLIC5"],
	"Secretory tubules and mesonephric remnant cells": ["GJA1", "CDH2", "AQP3", "CALB2"],
	"Female stromal cells": ["PDGFRA"],
	"Gonadal stromal cells": ["PDGFRA", "GATA4", "ARX"],
	"Mesonephric/fallopian tube stromal cells": ["PDGFRA", "GATA2"],
	"Proliferating stromal cells": ["PDGFRA", "MKI67"],
	"Mesonephric/epididymal/rete testis stromal cells": ["SULT1E1"],
	"Ovarian surface epithelium cells": ["WT1"],
	"Female mesonephros": ["ANXA2", "ANXA5"],
	"Ovary (f.MTs, f.WD, and f.MD)": ["ANXA1", "ANXA3"],
	"Female rete ovarii cells": ["PAX8", "PCP4", "BCAM", "SPRR2F"],
	"Rete ovarii epithelial cells": ["PAX8"],
	"f.MTs, f.WD, and f.MD in females": ["PAX8", "KRT19"],
    "Fetal Germ cells": ["POU5F1", "DDX4"]   
}

Taelman_meiosis = {
    "Epithelial rete testis": "somatic",
    "Sertoli": "somatic",
    "Endothelial cells": "somatic",
    "Smooth muscle cells": "somatic",
	"Immune cells": "somatic",
    "Epithelial cells from mesonephroi/epididymis": "somatic",
	"Stromal cells from mesonephroi/epididymis": "somatic",
	"Stromal cells from the testis": "somatic",
	"Leydig cells": "somatic",
	"Epithelial cells of the MTs and WD": "somatic",
	"Secretory tubules": "somatic",
	"Proliferating Sertoli cells, gonadal stromal cells, and mesonephric/epididymal stromal cells": "somatic",
	"Coelomic/ovarian surface epithelium": "somatic",
	"Pre-granulosa cells": "somatic",
	"Tubular cells of the f.MTs, f.MD, and f.WD": "somatic",
	"Mesonephric podocytes": "somatic",
	"Secretory tubules and mesonephric remnant cells": "somatic",
	"Female stromal cells": "somatic",
	"Gonadal stromal cells": "somatic",
	"Mesonephric/fallopian tube stromal cells": "somatic",
	"Proliferating stromal cells": "somatic",
	"Mesonephric/epididymal/rete testis stromal cells": "somatic",
	"Ovarian surface epithelium cells": "somatic",
	"Female mesonephros": "somatic",
	"Ovary (f.MTs, f.WD, and f.MD)": "somatic",
	"Female rete ovarii cells": "somatic",
	"Rete ovarii epithelial cells": "somatic",
	"f.MTs, f.WD, and f.MD in females": "somatic",
    "Fetal Germ cells": "meiotic"
}
"""
Irie = {
    "PGC_migra": ["TCL1A", "CXCR4"], 
    "PGC_mito": ["SUSD2", "RASSF2"],
    "PGC_mito_arrest": ["PIWIL1", "PIWIL2", "DNAJA4", "DDX4", "PIWIL4"]
}

Irie_meiosis = {
    "PGC_migra": "meiotic", 
    "PGC_mito": "meiotic",
    "PGC_mito_arrest": "meiotic"
}
"""
Irie = {
    "4i_ESCs": ["POU5F1", "DNMT3B"],
    "PGCLCs_N3+": ["SOX17", "PRDM1", "NANOS3", "BRDT", "BEND4", "KLF8", "HDAC4", "POU5F1"],
    "PGCLCs_DM+": ["CDH5", "DMRT1", "TCL1A", "SUSD2", "SOX17", "POU5F1", "TFCP2L1"],
    "PGCLCs_DZ+": ["DAZL", "PIWIL1", "PIWIL2", "DMRT1"]
}

Irie_meiosis = {
    "4i_ESCs": "meiotic",
    "PGCLCs_N3+": "meiotic",
    "PGCLCs_DM+": "meiotic",
    "PGCLCs_DZ+": "meiotic"
}

Seita = {
    "Cardiac lineage": ["CDKN1A", "EDN1","TGFB1","SPARC","SHC1",
    "CITED2","FOXF1","HEG1","MMP2","FN1","ECE1","PTN","GLI3",
    "ZFP36L1","PDLIM1","FKBP1A","FLRT2","COL5A1","PTCH1","CHD7",
    "FLNA","EDN1","EDNRA","HAND1","ECE1"],
    
    "Endoderm lineage": ["EOMES","SOX17","LHX1","NOG","DKK1","MIXL1",
    "VTN","COL4A2","MMP15","LAMB1","NODAL","MIXL1"],
    
    "Macrophage": ["CD84","NCF1","NCF2","CLEC10A","SLA","ADAR",
    "IFIT3","IFIH1","PYCARD","PSTPIP1","CLEC5A","TNFAIP8L2",
    "LYN","SPI1","TGFB1","TICAM1","MAPK14","MTDH","HCK","IRF3",
    "AKT1","CD14","TLR4","TLR2","MAPK3","PYCARD","FCER1G","SLC11A1",
    "IL2RB","BTK","CYBA","CD47","FCGR1A","DOCK2","IL2RG","MERTK",
    "DNM2"],
    
    "Endothelial lineage": ["ACVRL1","COL18A1","NRP1","ROBO4",
    "NRP2","FLT1","FLT4","TSPAN12","SOX18","KDR","VAV3","EGFL7",
    "UNC5B","CEMIP2","HSPG2","RHOB","COL13A1","SHC1","CTNND1",
    "ICAM2","THY1","F11R","PDLIM1","CDH5","NRP1","RAMP2","EGFL7",
    "TGFB1","FZD4","TIE1","CAV1","QKI","RASIP1"],
    
    "Apoptotic": ["RARG","SH3KBP1","UBE2D3","PDCD5","PTEN",
    "SRA1","UBE2Z","BBC3","ING4","BCL7C","RACK1","KIF1B","JAK2",
    "RHOB","DDB1","ALDH1A3","RABEP1","TRAF4","AX1BP1","ELMO2",
    "BAX","TAF6","BCL2L1","PDIA3","GSK3B","TGFB2","TGFB1","PDPK1",
    "PARP2","RB1CC1","BAX","INHBA","JAK2","VCP","MAGED1","SRA1",
    "ETS1","MBTD1","ING2","EPC1","EPC2","JAK2","CCNL1","PAK2",
    "ZBTB7A","RBM5","EGR1","NME2","FBXO11","NME1","MORF4L1",
    "TRAF4","GAS1","SGK1","BCL2L1"],
    
    "iPSCs": ["ANAPC13","CDCA3","CCAR1","CKS1B","PCNP","PPP1CC",
    "CHAF1A","MIS18BP1","CDC26","NUF2","TLK1","HMG20B","OIP5",
    "CASP8AP2","CDC26","NUF2","OIP5","HELLS","CENPV","NUDC",
    "TIPIN","CENPW","DYNLT1","BOD1","SLF1","NSMCE4A","MYC","NSMCE1"],
    
    "PGCLCs_d2": ["EOMES","TWSG1","BMPR2","HAND1","NF2","BMP7",
    "WNT3","WLS","USP14","KLHDC3","VCP","FBXW5","ZFAND2A","HUWE1",
    "UBE2A","AXIN2","RAD23B","PSMB10","PSMA7","AURKA","CD2AP","PCNP",
    "DDB1","CDC34","HECTD1","PSMC3","PSMD2","TRIM2","PCBP2"],
    
    "PGCLCs_d6": ["KDM5B","KDM3B","KMT2A","CHD9","PHF20","ABRAXAS1",
    "KMT2C","CHD6","JMJD1C","ARID4B","NR3C1","PHF8","NSD3","BABAM2",
    "HMGN4","SOX15","SMYD2","HMG20A","SUDS3","BRD1","GPX4","USP3",
    "SIRT1","RSBN1","KANSL1","DPPA3","HAT1","RCOR1","ATF7IP",
    "KMT2A","MGMT","ATRX","BAZ2A","BMI1","KLF10","RIF1","VPS72",
    "STAT3","KIT","NANOG","RBPJ","SOX4"]
    
}

Seita_meiosis = {
    "Cardiac lineage": "somatic",
    "Endoderm lineage": "somatic",
    "Macrophage": "somatic",
    "Endothelial lineage": "somatic",
    "Apoptotic": "somatic",
    "iPSCs": "meiotic",
    "PGCLCs_d2": "meiotic",
    "PGCLCs_d6": "meiotic"  
}

Huang = {
    "Immature Sertoli cell": ["AMH","INHA","CLU"],
	"Mature Sertoli cell": ["FATE1", "HMGN5", "CLDN11"],
	"Peritubular myoid cell": ["DCN", "PDGFRA", "PTCH1"],
	"Leydig cell": ["STAR", "HSD3B1", "INSL3"],
	"Macrophage": ["CSF1R", "CD74", "MAFB"],
	"Endothelial cells": ["PECAM1", "CLDN5", "ECSCR"], 
	"Natural killer cell": ["PTPRC", "NKG7", "KLRF1", "CD94"],
	"SPG_undiff": ["SALL4", "ZBTB16", "ELAVL2"],
	"SPG": ["NR6A1", "FGFR3", "FMR1"],
	"SPC_prelep_zyg": ["TEX12", "SMCHD1", "PRSS50"],
	"SPC_pach_secondary": ["SYCE2", "BOLL", "MLH1"],
	"STD_r": ["DNAH14", "IZUMO4", "NKAPL"],
	"STD_e": ["TEX29", "FBXO24", "ACRV1"],
	"Sperm": ["ODF1", "ODF2", "OAZ3", "AKAP4"]
}

Huang_meiosis = {
    "Immature Sertoli cell": "somatic",
	"Mature Sertoli cell": "somatic",
	"Peritubular myoid cell": "somatic",
	"Leydig cell": "somatic",
	"Macrophage": "somatic",
	"Endothelial cells": "somatic", 
	"Natural killer cell": "somatic",
	"SPG_undiff": "meiotic",
	"SPG": "meiotic",
	"SPC_prelep_zyg": "meiotic",
	"SPC_pach_secondary": "meiotic",
	"STD_r": "meiotic",
	"STD_e": "meiotic",
	"Sperm": "meiotic"
}

Sosa = {
    "Mesoderm": ["TBX4", "HAND1", "DPPA3", "POU5F1", "KDR", "TAL1"],
    "Trophoblast / amnion": ["TFAP2C", "HAND1"],
    "Stem-like cells": ["POU5F1", "NANOG", "SOX2"],
    "PGCLCs": ["POU5F1", "NANOG", "SOX2", "DPPA3", "TFAP2C", "PRDM1", "PRDM14", "DAZL", "DDX4", "NANOS3"]
}

Sosa_meiosis = {
    "Mesoderm": "somatic",
    "Trophoblast / amnion": "somatic",
    "Stem-like cells": "somatic",
    "PGCLCs": "meiotic"
}

Merrick = {
    "Germ cells": ["DAZL", "DDX4", "FGFR3", "UTF1"],
    "Mesothelial": ["UPK3B"],
    "Gonadal somatic": ["GATA4", "LHX9", "NR5A1"],
    "Supporting": ["WNT6"],
    "pre-Granulosa": ["IRX3", "FOXL2"],
    "Sertoli": ["AMH", "SOX9"],
    "Gonadal interstitial": ["ARX", "TCF21"],
    "Mesenchymal": ["PDGFRA", "DCN"],
    "Extragonadal": ["GATA2", "NR2F1"],
    "PV": ["PDGFRB"],
    "SMC": ["MYH11"],
    "Immune": ["PTPRC"],
    "Epithelial": ["PAX8", "EPCAM"],
    "Erythroid": ["HBA1"],
    "Neural": ["ASCL1"],
    "Macrophages": ["CD14", "CD163", "S100A4", "TYROBP", "LYZ", "RGS1"],
    "Endothelial cells": ["CDH5", "VWF", "EPAS1", "TGFBR2", "NOSTRIN", "PALMD", "POSTN"],
    "Myoid cells": ["ACTA2", "MYH11", "TPM1", "TPM2", "TPM4", "MYL9"],
    "Leydig cells": ["DLK1", "IGF1", "IGF2", "IGFBP5", "HSD17B3", "CFD", "INSL3", "TCF21", "ARX"],
    "Somatic cells": ["VIM"],
    "SSC_state0": ["PIWIL4", "PHGDH", "SLC25A22", "ICA1L", "PPP1R36", "MAGEB1", "EGR4", "MSL3", "TSPAN33"],
    "SSC": ["UTF1", "ID4", "FGFR3", "GFRA1", "ETV5", "ZBTB16", "MAGEA4", "TCF3", "L1TD1"],
    "SPG": ["TBX3", "HOXA3"],
    "SPG_diff": ["KIT", "DMRT1", "MKI67", "SOHLH1", "SOHLH2", "STRA8", "MAGEA4"],
    "SCT_prelep": ["STRA8", "REC8", "DMRT1", "SOX4"],
    "SCT_early": ["CHEK1", "BRCA1", "SPO11", "DMC1", "ATM", "SYCP1", "SYCP2", "SYCP3", "MAGEA4"],
    "SCT_pach_dyp": ["SYCP3", "MLH3", "MAGEA4", "CCNA1"],
    "STD_r": ["SPAG6"],
    "STD_e": ["ZPBP", "ZPBP2", "CAMK4", "DNAH6", "DNAH7", "DNAH14", "CATSPER1", "CATSPER4", "CREM", "MYO1D"],
    "Sperm": ["TNP1", "TNP2", "PRM2", "HOOK1", "SPATA7", "SPATA32", "SPATA33", "SPATA12", "SPATA18", "SPATA20", "PRM3"]  
}

Merrick_meiosis = {
    "Germ cells": "meiotic",
    "Mesothelial": "somatic",
    "Gonadal somatic": "somatic",
    "Supporting": "somatic",
    "pre-Granulosa": "somatic",
    "Sertoli": "somatic",
    "Gonadal interstitial": "somatic",
    "Mesenchymal": "somatic",
    "Extragonadal": "somatic",
    "PV": "somatic",
    "SMC": "somatic",
    "Immune": "somatic",
    "Epithelial": "somatic",
    "Erythroid": "somatic",
    "Neural": "somatic",
    "Macrophages": "somatic",
    "Endothelial cells": "somatic",
    "Myoid cells": "somatic",
    "Leydig cells": "somatic",
    "Somatic cells": "somatic",
    "SSC_state0": "meiotic",
    "SSC": "meiotic",
    "SPG": "meiotic",
    "SPG_diff": "meiotic",
    "SCT_prelep": "meiotic",
    "SCT_early": "meiotic",
    "SCT_pach_dyp": "meiotic",
    "STD_r": "meiotic",
    "STD_e": "meiotic",
    "Sperm": "meiotic"  
}

oocytes = {
    "STD": ["ZPBP", "ZPBP2", "CAMK4", "DNAH6", "DNAH7", "DNAH14", "CATSPER1", "CATSPER4", "CREM", "MYO1D"]
}

oocytes_meiosis = {
    "STD": "meiotic"
}

embryo_female = {
    "PGCs mitotic": ["POU5F1", "SOX2", "UTF1", "DAZL", "DDX4", "CENPF", "SALL4"],
    "PGCs early meiotic": ["POU5F1", "SOX2", "UTF1", "DAZL", "DDX4", "CENPF", "STRA8", "SYCP1", "SYCP3", "REC8", "SMC1B", "DUSP9", "WDR89", "DPPA5A"],
    "PGCs late meiotic": ["DAZL", "DDX4", "CENPF", "STRA8", "SYCP1", "SYCP3", "TEX14", "MAEL", "TEX101", "TAF7L"],
    "Pregranulosa": ["WNT4", "WNT6", "KCTD14"],
    "Mesothelial": ["LHX9", "UPK3B", "ALDH1A2"],
    "Endothelial": ["PECAM1", "FLT1", "KDR"],
    "Interstitial": ["COL1A2", "COL1A1", "BGN"],
    "Immune cells": ["CD52", "CAR2"],
    "Erythroid cells": ["ALAS2", "ALAD"]
}

embryo_female_meiosis = {
    "PGCs mitotic": "meiotic",
    "PGCs early meiotic": "meiotic",
    "PGCs late meiotic": "meiotic",
    "Pregranulosa": "somatic",
    "Mesothelial": "somatic",
    "Endothelial": "somatic",
    "Interstitial": "somatic",
    "Immune cells": "somatic",
    "Erythroid cells": "somatic"
}

comb_list = {
    "SPG": "SPG",
    "SCT_lep": "SCT",
    "SCT_zyg": "SCT",
    "SCT_pach": "SCT",
    "SCT_dyp": "SCT",
    "STD_early": "STD",
    "STD_late": "STD",
    "PMCs": "PMCs",
    "Fibrotic_PMCs": "Fibrotic_PMCs",
    "Sertoli_cells": "Sertoli_cells",
    "Leydig_cells": "Leydig_cells",
    "Endothelial_cells": "Endothelial_cells",
    "Macrophages": "Macrophages",
    "SSC": "SPG",
    "SPG_prog": "SPG",
    "SPG_early_diff": "SPG",
    "SPG_late_diff": "SPG",
    "SCT_prelep": "SCT",
    "SCT_lep_zyg": "SCT",
    "SCT_pach": "SCT",
    "SCT_dyp": "SCT",
    "STD_early_r": "STD",
    "STD_mid_r": "STD",
    "STD_late_r": "STD",
    "PMCs": "PMCs",
    "Fibrotic_PMCs": "Fibrotic_PMCs",
    "Sertoli_cells": "Sertoli_cells",
    "Leydig_cells": "Leydig_cells",
    "Endothelial_cells": "Endothelial_cells",
    "Macrophages": "Macrophages",
    "Sertoli_cells": "Sertoli_cells",
    "PSMCs": "PSMCs",
    "Endothelial cells": "Endothelial cells",
    "Macrophages": "Macrophages",
    "Leydig_cells": "Leydig_cells",
    "SPG_undiff": "SPG",
    "SPG_diff": "SPG",
    "SCT_lep": "SCT",
    "SCT_zyg": "SCT",
    "SCT_pach": "SCT",
    "STD_early_r": "STD",
    "STD_late_r": "STD",
    "STD_e": "STD",
    "Epithelial rete testis": "Epithelial rete testis",
    "Sertoli": "Sertoli",
    "Endothelial cells": "Endothelial cells",
    "Smooth muscle cells": "Smooth muscle cells",
	"Immune cells": "Immune cells",
    "Epithelial cells from mesonephroi/epididymis": "Epithelial cells from mesonephroi/epididymis",
	"Stromal cells from mesonephroi/epididymis": "Stromal cells from mesonephroi/epididymis",
	"Stromal cells from the testis": "Stromal cells from the testis",
	"Leydig cells": "Leydig cells",
	"Epithelial cells of the MTs and WD": "Epithelial cells of the MTs and WD",
	"Secretory tubules": "Secretory tubules",
	"Proliferating Sertoli cells, gonadal stromal cells, and mesonephric/epididymal stromal cells": "Proliferating Sertoli cells, gonadal stromal cells, and mesonephric/epididymal stromal cells",
	"Coelomic/ovarian surface epithelium": "Coelomic/ovarian surface epithelium",
	"Pre-granulosa cells": "Pre-granulosa cells",
	"Tubular cells of the f.MTs, f.MD, and f.WD": "Tubular cells of the f.MTs, f.MD, and f.WD",
	"Mesonephric podocytes": "Mesonephric podocytes",
	"Secretory tubules and mesonephric remnant cells": "Secretory tubules and mesonephric remnant cells",
	"Female stromal cells": "Female stromal cells",
	"Gonadal stromal cells": "Gonadal stromal cells",
	"Mesonephric/fallopian tube stromal cells": "Mesonephric/fallopian tube stromal cells",
	"Proliferating stromal cells": "Proliferating stromal cells",
	"Mesonephric/epididymal/rete testis stromal cells": "Mesonephric/epididymal/rete testis stromal cells",
	"Ovarian surface epithelium cells": "Ovarian surface epithelium cells",
	"Female mesonephros": "Female mesonephros",
	"Ovary (f.MTs, f.WD, and f.MD)": "Ovary (f.MTs, f.WD, and f.MD)",
	"Female rete ovarii cells": "Female rete ovarii cells",
	"Rete ovarii epithelial cells": "Rete ovarii epithelial cells",
	"f.MTs, f.WD, and f.MD in females": "f.MTs, f.WD, and f.MD in females",
    "Fetal Germ cells": "Fetal Germ cells",
    "4i_ESCs": "ESCs",
    "PGCLCs_N3+": "PGCLCs",
    "PGCLCs_DM+": "PGCLCs",
    "PGCLCs_DZ+": "PGCLCs",
    "Cardiac lineage": "Cardiac lineage",
    "Endoderm lineage": "Endoderm lineage",
    "Macrophage": "Macrophage",
    "Endothelial lineage": "Endothelial lineage",
    "Apoptotic": "Apoptotic",
    "iPSCs": "iPSCs",
    "PGCLCs_d2": "PGCLCs",
    "PGCLCs_d6": "PGCLCs",
    "Immature Sertoli cell": "Immature Sertoli cell",
	"Mature Sertoli cell": "Mature Sertoli cell",
	"Peritubular myoid cell": "Peritubular myoid cell",
	"Leydig cell": "Leydig cell",
	"Macrophage": "Macrophage",
	"Endothelial cells": "Endothelial cells", 
	"Natural killer cell": "Natural killer cell",
	"SPG_undiff": "SPG",
	"SPG": "SPG",
	"SPC_prelep_zyg": "SCT",
	"SPC_pach_secondary": "SCT",
	"STD_r": "STD",
	"STD_e": "STD",
	"Sperm": "STD",
    "Mesoderm": "Mesoderm",
    "Trophoblast / amnion": "Trophoblast / amnion",
    "Stem-like cells": "Stem-like cells",
    "PGCLCs": "PGCLCs",
    # Merrick
    "Germ cells": "Fetal Germ cells",
    "Mesothelial": "Mesothelial",
    "Gonadal somatic": "Gonadal somatic",
    "Supporting": "Supporting",
    "pre-Granulosa": "pre-Granulosa",
    "Sertoli": "Sertoli",
    "Gonadal interstitial": "Gonadal interstitial",
    "Mesenchymal": "Mesenchymal",
    "Extragonadal": "Extragonadal",
    "PV": "PV",
    "SMC": "SMC",
    "Immune": "Immune",
    "Epithelial": "Epithelial",
    "Erythroid": "Erythroid",
    "Neural": "Neural",
    "Macrophages": "Macrophages",
    "Endothelial cells": "Endothelial cells",
    "Myoid cells": "Myoid cells",
    "Leydig cells": "Leydig cells",
    "Somatic cells": "Somatic cells",
    "SSC_state0": "SPG",
    "SSC": "SPG",
    "SPG": "SPG",
    "SPG_diff": "SPG",
    "SCT_prelep": "SCT",
    "SCT_early": "SCT",
    "SCT_pach_dyp": "SCT",
    "STD_r": "STD",
    "STD_e": "STD",
    "Sperm": "STD",
    # oocytes
    "STD": "STD",
    # embryo_female
    "PGCs mitotic": "PGCs",
    "PGCs early meiotic": "PGCs",
    "PGCs late meiotic": "PGCs",
    "Pregranulosa": "Pregranulosa",
    "Mesothelial": "Mesothelial",
    "Endothelial": "Endothelial",
    "Interstitial": "Interstitial",
    "Immune cells": "Immune cells",
    "Erythroid cells": "Erythroid cells"
}
