---
layout: post
date: 2016-01-09
title: "Why you should use ete for tree exploration and visualisation in python"
tags:
    - python
    - notebook
    - phylogenetics
---

If you work with trees (be it phylogenetics or not) and you regularly use python, you have probably used or heard about one of the following packages: [Bio.phylo](https://github.com/biopython/biopython/tree/master/Bio/Phylo), [dendropy](https://pythonhosted.org/DendroPy/) and [ete](etetoolkit.org). 

While each one of those packages has its own unique strengths and weaknesses., I particulary like the **ETE** module. Here is why !

<!--more-->

This post is based on one of my past presentation at [monbug](http://www.monbug.ca/). The github repository with all the file can be found here : https://github.com/maclandrol/monbug_ete . You could use [nbviewer](http://nbviewer.ipython.org/github/maclandrol/monbug_ete/blob/master/mon_bug_sep3.ipynb) to view the notebook if you prefer.

### What's ETE ??
ETE is a python **E**nvironment for **T**ree **E**xploration created by _Jaime Huerta-Cepas_.

It's a framework that assists in the manipulation of any type of hierarchical tree (this include reading, writing, visualisation, annotation, etc). The current version is ```ete3```

### Installation 

You can install ETE simply by using pip : ```pip install ete3```. Check this link for more detail about optionnal/unmet dependencies : http://etetoolkit.org/download/


### Quick introduction to the API 

A great in-depth tutorial for working with tree data structure in ete is provided by the author here : http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html 

I'll be making a light introduction to the API here, but I really recommend you to read the official doc! Let's take a quick glance at the available tree data structure in ete : 

**In [58]:**

{% highlight python %}
import ete3
import inspect
print([x[0] for x in inspect.getmembers(ete3, inspect.isclass) if x[0].endswith('Tree')])

{% endhighlight %}

    ['ClusterTree', 'EvolTree', 'NexmlTree', 'PhyloTree', 'PhyloxmlTree', 'Tree']

As you can see, you have a basic tree data structure (```Tree```) and more specific tree structures like ```PhyloTree``` for __phylogenetics__

#### ETE can read tree from a str or a file 

**In [59]:**

{% highlight python %}
from ete3 import Tree

rand_newick = "((((a,b), c), d), (e,f));"
rand_tree = "rand_tree"
with open(rand_tree, 'w') as TREE:
    TREE.write(rand_newick)

# Reading tree t1 and t2
# FILL HERE
t1 = Tree(rand_newick)
t2 = Tree(rand_tree)

{% endhighlight %}
 

#### A tree is a Node ===> the root is a Node, so are all its descendents 

**In [61]:**

{% highlight python %}
# Print t1 then t2 in ASCII
print(t1.get_ascii())
print(t2)
{% endhighlight %}

    
                /-a
             /-|
          /-|   \-b
         |  |
       /-|   \-c
      |  |
    --|   \-d
      |
      |   /-e
       \-|
          \-f
    
                /-a
             /-|
          /-|   \-b
         |  |
       /-|   \-c
      |  |
    --|   \-d
      |
      |   /-e
       \-|
          \-f

 
### You can add features to nodes 

**In [62]:**

{% highlight python %}
# Traverse tree t1 and add a feature at internal_node
from numpy import random

# Traverse : levelorder, preorder, postorder
for node in t1.traverse("levelorder"):
    if node.is_leaf():
        # add a features : randomness
        node_rand = random.randint(10)
        # FILL CODE
        node.add_features(randomness=node_rand)
{% endhighlight %}
 
### Features are just attributes
 

**In [63]:**

{% highlight python %}
# print t1 again with features : name and randomness
print(t1.get_ascii(attributes=['name', 'randomness']))

# Iterate over all the leaves
# t1.iter_leaves()
print("\n")
for node in t1:
    print("%s\t%s"%(node.name, node.randomness))

{% endhighlight %}

    
                /-a, 8
             /-|
          /-|   \-b, 1
         |  |
       /-|   \-c, 9
      |  |
    --|   \-d, 3
      |
      |   /-e, 9
       \-|
          \-f, 3
    
    
    a	8
    b	1
    c	9
    d	3
    e	9
    f	3

 
### You can search by features 

**In [64]:**

{% highlight python %}
# search by features
print(t1.search_nodes(randomness=8))
print(t1.search_nodes(name='a'))

# shortcut for node name :)
t1&'a'
{% endhighlight %}

    [Tree node 'a' (-0x7ffff810443aa570)]
    [Tree node 'a' (-0x7ffff810443aa570)]

    Tree node 'a' (-0x7ffff810443aa570)


 
### A quick list of useful methods 

**In [65]:**

{% highlight python %}
print(t1)
# get sister node =====> get_sisters()
sister = (t1&'a').get_sisters()
print("\nSISTER a : ")
print(sister)

{% endhighlight %}

    
                /-a
             /-|
          /-|   \-b
         |  |
       /-|   \-c
      |  |
    --|   \-d
      |
      |   /-e
       \-|
          \-f
    
    SISTER a : 
    [Tree node 'b' (0x7efbbc55ab0)]


**In [66]:**

{% highlight python %}
print(t1)
# get children  ===< get_children()
root_children = t1.get_children()
print("\n\nFIRST CHILD OF ROOT")
print(root_children[0])
{% endhighlight %}

    
                /-a
             /-|
          /-|   \-b
         |  |
       /-|   \-c
      |  |
    --|   \-d
      |
      |   /-e
       \-|
          \-f
    
    
    FIRST CHILD OF ROOT
    
             /-a
          /-|
       /-|   \-b
      |  |
    --|   \-c
      |
       \-d


**In [67]:**

{% highlight python %}
print(t1)
# Get LCA of multiple node ===> get_common_ancestor()
lca = t1.get_common_ancestor(['a', 'b'])
print("\n\nLCA (a, b) : ")
print(lca)
{% endhighlight %}

    
                /-a
             /-|
          /-|   \-b
         |  |
       /-|   \-c
      |  |
    --|   \-d
      |
      |   /-e
       \-|
          \-f
    
    
    LCA (a, b) : 
    
       /-a
    --|
       \-b


**In [68]:**

{% highlight python %}
# RF distance
# t1 and t2 were the same tree ... 
rf = t1.robinson_foulds(t2)
print("\n\nRF DISTANCE between t1 and t2 :")
print(rf[0])
{% endhighlight %}

    
    
    RF DISTANCE between t1 and t2 :
    0


**In [69]:**

{% highlight python %}
# prune tree to list of leaf: ===> prune
leaf_to_keep = ['a', 'c', 'd', 'f']
print("\n\nPRUNE t1 to [%s]"%", ".join(leaf_to_keep))
print("\n\nBEFORE PRUNING : ")
print(t1)

# FILL CODE
t1.prune(leaf_to_keep)

print("\n\nAFTER PRUNING : ")
print(t1)
{% endhighlight %}

    
    
    PRUNE t1 to [a, c, d, f]
    
    
    BEFORE PRUNING : 
    
                /-a
             /-|
          /-|   \-b
         |  |
       /-|   \-c
      |  |
    --|   \-d
      |
      |   /-e
       \-|
          \-f
    
    
    AFTER PRUNING : 
    
             /-c
          /-|
       /-|   \-a
      |  |
    --|   \-d
      |
       \-f

 
## PART 2 : Visualization 
 
### Example 1 : simple tree visualisation

Data : a random tree with random branches

- Tree rendering
- Tree Style
 

**In [71]:**

{% highlight python %}
from ete3 import Tree

# Generate a random tree (yule process) with leaf name defined by leave_names
t = Tree()
t.populate(8, names_library=list('ABCDEFGHIJKL'), random_branches=True)

print(t.get_ascii(attributes=['name', 'support'], show_internal=True))

{% endhighlight %}

    
                              /-G, 0.4793609239053189
         /, 0.11319458200540156
        |                    |                   /-F, 0.5340382546548448
        |                     \, 0.520949167241861
    -, 1.0                                       \-E, 0.8982207041579878
        |
        |                     /-L, 0.2768295431285249
         \, 0.32620620577540727
                             |                     /-K, 0.5017351222348678
                              \, 0.07320524397375128
                                                  |                   /-J, 0.1420801709839633
                                                   \, 0.931415565866013
                                                                     |                   /-I, 0.05555721516881551
                                                                      \, 0.875123316663256
                                                                                         \-H, 0.8108877272041851

 
### Show tree in a gui
 

**In [72]:**

{% highlight python %}

t.show()
{% endhighlight %}
 
### Render tree (supported format : png, pdf and svg)
### Notebook support !!!!
 

**In [74]:**

{% highlight python %}
from IPython.display import Image

# render tree (supported format : png, pdf and svg)
#t.render('tree.png', dpi=200)
#Image('tree.png')

# render and show directly in the notebook
# FILL CODE
t.render('%%inline', dpi=200)
{% endhighlight %}



 
![png]({{ site.baseurl }}/public/images/mon_bug_oct_files/mon_bug_oct_29_0.png) 


 
### ==> This is what make ETE different 

**In [75]:**

{% highlight python %}
# TreeStyle to change the display options

from ete3 import TreeStyle

ts = TreeStyle()

ts.show_branch_length = True # show branch length
ts.show_branch_support = True # show support

# rotate  tree
ts.rotation = -30
t.render('%%inline', tree_style=ts)
{% endhighlight %}



 
![png]({{ site.baseurl }}/public/images/mon_bug_oct_files/mon_bug_oct_31_0.png) 



**In [76]:**

{% highlight python %}
# Draw a circular tree

ts.rotation = 0
ts.mode = "c" #circular mode (default = 'r' for rectangular ?)
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 180
t.render('%%inline', tree_style=ts, w=500)
{% endhighlight %}



 
![png]({{ site.baseurl }}/public/images/mon_bug_oct_files/mon_bug_oct_32_0.png) 


 
### The wonder of  "faces'' 

**In [77]:**

{% highlight python %}
# Adding face to Tree
from ete3 import faces
print([f for f in dir(faces) if 'Face' in f])
Image('face_positions.png')
{% endhighlight %}

    ['AttrFace', 'BarChartFace', 'CircleFace', 'DynamicItemFace', 'Face', 'ImgFace', 'OLD_SequenceFace', 'PieChartFace', 'ProfileFace', 'RandomFace', 'RectFace', 'SeqMotifFace', 'SequenceFace', 'SequencePlotFace', 'StackedBarFace', 'StaticItemFace', 'TextFace', 'TreeFace']




 
![png]({{ site.baseurl }}/public/images/mon_bug_oct_files/mon_bug_oct_34_1.png) 



**In [78]:**

{% highlight python %}
# Adding a Circleface to our current tree

def layout(node):
    if node.is_leaf():
     
        # Create a sphere face with a random size
        C = faces.CircleFace(radius=random.randint(5, 50), color="RoyalBlue", style="sphere")
        
        # Make it transparent 
        C.opacity = 0.3
       
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")

# link out layout to the style
ts.layout_fn = layout

#show
t.render('%%inline', tree_style=ts, w=700)
{% endhighlight %}



 
![png]({{ site.baseurl }}/public/images/mon_bug_oct_files/mon_bug_oct_35_0.png) 



**In [79]:**

{% highlight python %}
# With Face, You can actually make things like this : (treeception)
Image('tree_faces.png')
{% endhighlight %}



 
![png]({{ site.baseurl }}/public/images/mon_bug_oct_files/mon_bug_oct_36_0.png) 


 
### Example 2 : Duplication|Loss history of a gene familly

Data : genetree newick where I already add a feature (**state**), where :

- states = 1 ==> internal node with duplication
- states = 0 ==> speciation node
 

**In [80]:**

{% highlight python %}
from ete3 import Tree
t = Tree('annoted_trees', format=2)
print(t.get_ascii(show_internal=True, attributes=['name', 'states']))
{% endhighlight %}

    
          /-Dre_1, 0
       /, 0
      |  |   /-Cfa_1, 0
      |   \, 0
    -, 1     \-Hsa_1, 0
      |
      |   /-Dre_2, 0
       \, 0
          \-Cfa_2, 0


**In [81]:**

{% highlight python %}
from ete3 import Tree, faces, TreeStyle
import utils

# Creates a layout function
def mylayout(node):
    # If node is a leaf, its scientific name
    if node.is_leaf():
       
        # add a face for scientific name
        longNameFace = faces.TextFace(utils.get_scientific_name(node))
        faces.add_face_to_node(longNameFace, node, column=1)

        # add an image Face
        node.img_style["size"] = 0
        image = utils.get_image(node.name)
        faces.add_face_to_node(faces.ImgFace(image), node, column=0, aligned=True)
    
    #If node is an internal node
    elif int(node.states) == 1:
        # Sets the style of internal nodes
        node.img_style["size"] = 6
        node.img_style["shape"] = "square"
        node.img_style["fgcolor"] = "green"
    else :
        # Sets the style of internal nodes
        node.img_style["size"] = 6
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "red"
        

# And, finally, display the tree using my own layout function
ts = TreeStyle()
ts.show_leaf_name = False
ts.layout_fn = mylayout

t.render("%%inline", dpi=600, tree_style = ts)


{% endhighlight %}



 
![png]({{ site.baseurl }}/public/images/mon_bug_oct_files/mon_bug_oct_39_0.png) 


 
### Example 3 : Phylogenetic tree, sequence and information content

Data :
    - An alignment
    - A tree constructed using that alignment
    (Actually those two were randomly generated)
 

**In [82]:**

{% highlight python %}
from ete3 import PhyloNode, SequenceFace, faces, TreeStyle
from Bio import AlignIO
from Bio import Alphabet
from Bio.Align import AlignInfo
from utils import show_file

alignement = "alignement.fasta"
tree = "phylotree.nw"

# Open tree and link alignement to it
t = PhyloNode(tree)
t.link_to_alignment(alignement)
show_file(alignement)
show_file(tree)
{% endhighlight %}

    >A
    MAEIPDETIQQFMALT---SNIAVQYLSEFGDLNEALNSY
    >B
    MAEIPDATIQQFMALTNVSHNIAVQY--EFGDLNEALNSY
    >C
    MAEIPDATIQ----LTNVSHNIAVQYLSEFGDLNEALNSY
    >D
    MAEAPDETIQQFMALTNVSHNIAVQYLSEFGDLNEAL---
    (A,(D,(B,C)));


**In [83]:**

{% highlight python %}
# Compute Information content with Biopython
align = AlignIO.read(alignement, 'fasta', alphabet=Alphabet.Gapped(Alphabet.IUPAC.protein))
summary_info = AlignInfo.SummaryInfo(align)        
total_ic_content = summary_info.information_content()
ic_content = summary_info.ic_vector.values()

# Set TreeStyle
ts = TreeStyle()
ts.branch_vertical_margin = 10
ts.allow_face_overlap = False
ts.show_scale = False
ts.show_leaf_name = False

# Align ic plot to TreeStyle header
ic_plot = faces.SequencePlotFace(ic_content, fsize=10, col_width=14, header="Information Content", kind='bar', ylabel="ic")
ts.aligned_header.add_face(ic_plot, 1) 

#t.add_face(ic_plot,1)
t.render("%%inline", tree_style=ts, dpi=300)
{% endhighlight %}



 
![png]({{ site.baseurl }}/public/images/mon_bug_oct_files/mon_bug_oct_42_0.png) 


 
## THANKS

- you'll find the notebook here : http://nbviewer.ipython.org/github/maclandrol/
monbug_ete/blob/master/mon_bug_sep3.ipynb

- and the associated github repo here : https://github.com/maclandrol/monbug_ete

- ETE toolkit : http://etetoolkit.org/

- If you want to contribute (Cepas github repo): https://github.com/jhcepas/ete 

**In [None]:**

{% highlight python %}

{% endhighlight %}
