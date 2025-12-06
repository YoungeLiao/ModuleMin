# %%
from ete3 import Tree, TreeStyle, NodeStyle

# %%
def visualize_tree(tree_file):
    # Load the tree
    tree = Tree(tree_file, format=1)

    # Define a tree style
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True

    # Define a node style for clades
    def set_clade_style(node, color):
        nstyle = NodeStyle()
        nstyle["fgcolor"] = color
        nstyle["size"] = 10
        node.set_style(nstyle)

    # Identify and style clades
    for node in tree.traverse():
        if node.is_leaf():
            continue
        if len(node) > 1:  # Example condition to identify clades
            set_clade_style(node, "red")

    # Render the tree
    tree.show(tree_style=ts)

if __name__ == "__main__":
    tree_file = "../results/output_ref.tree"
    visualize_tree(tree_file)
# %%
