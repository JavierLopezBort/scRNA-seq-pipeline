meiotic_cells = ["germ cell", "oogonial cell", "epithelial cell", "pancreatic acinar cell", "thyroid follicular cell", "skeletal muscle satellite stem cell", "classical monocyte", "hepatocyte", "placental villous trophoblast", "monocyte", "endocrine cell"]
somatic_cells = ["macrophage", "myofibroblast cell", "endothelial cell", "Sertoli cell", "fibroblast", "blood vessel endothelial cell", "granulosa cell", "mesenchymal cell", "myoepithelial cell of mammary gland", "neuron", "pericyte", "supporting cell", "vascular associated smooth muscle cell", "endothelial cell of lymphatic vessel"]
comb_list = {
    "germ cell": "SPG", # PGCLCs # Fetal_germ_cells
    "oogonial cell": "SCT", #
    "epithelial cell": "STD", #
    "pancreatic acinar cell": "STD", # 
    "thyroid follicular cell": "STD", #
    "skeletal muscle satellite stem cell": "STD",
    "classical monocyte": "STD", #
    "hepatocyte": "PGCLCs", # somatic # PGCLCs
    "placental villous trophoblast": "PGCLCs",
    "monocyte": "STD",
    "endocrine cell": "STD",
    "macrophage": "macrophage",
    "myofibroblast cell": "myofibroblast cell",
    "endothelial cell": "endothelial cell",
    "Sertoli cell": "Sertoli cell",
    "fibroblast": "fibroblast",
    "blood vessel endothelial cell": "blood vessel endothelial cell",
    "granulosa cell": "granulosa cell",
    "mesenchymal cell": "mesenchymal cell",
    "myoepithelial cell of mammary gland": "myoepithelial cell of mammary gland",
    "neuron": "neuron",
    "pericyte": "pericyte",
    "supporting cell": "supporting cell",
    "vascular associated smooth muscle cell": "vascular associated smooth muscle cell",
    "endothelial cell of lymphatic vessel": "endothelial cell of lymphatic vessel"
}