# A simple code to pull NCBI taxonomy tree and save it in png format
# Based on the https://www.biostars.org/p/153788/#154006

from ete3 import NCBITaxa, AttrFace, TreeStyle
# This may take a while if NCBI taxonomy has not been loaded to the local system before
ncbi = NCBITaxa()
# A list of taxid of organisms
tree =
ncbi.get_topology([240159,215358,210632,56716,8187,8167,8168,80972,161767,80966,144197,8
175,244447,28829,8022,229290,69293,205130])
# custom layout: adds "rank" on top of branches, and sci_name as tip names
def my_layout(node):
if getattr(node, "rank", None):
rank_face = AttrFace("rank", fsize=7, fgcolor="indianred")
node.add_face(rank_face, column=0, position="branch-top")
if node.is_leaf():
sciname_face = AttrFace("sci_name", fsize=9, fgcolor="steelblue")
node.add_face(sciname_face, column=0, position="branch-right")
ts = TreeStyle()
ts.layout_fn = my_layout
ts.show_leaf_name = False
ts.show_scale = False
tree.render("tree2.png", tree_style=ts, w=1000, units="px")
