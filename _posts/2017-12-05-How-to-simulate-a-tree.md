---
layout: post
date: 2017-12-05
title: "How to simulate a phylogenetic tree ? (part 1)"
tags: ['phylogenetic', 'tree', 'bioinfo', 'python']
comments: true
math: true
---

Benchmarking phylogenetic algorithms often required the availability of a dataset of true evolutionary histories. As these true histories are hardly known, we rely instead on simulated datasets. In this series of post, I will present simple algorithms for the simulation of phylogenetic trees depicting the evolutionary history of a gene family. This first part will focus on the simulation of `species tree` under either a `pure birth` or a `birth-death` process.

<!--more-->

Simulations of phylogenetic trees are often conducted under a constant rate of evolution (molecular clock). The molecular clock can later be relaxed to allow variation in rates accross branches. In fact, the algorithms presented here will only consider elapsed time. Therefore, branch lengths will only represent time (no substitution rate). 

## Ultrametric binary tree with a pure birth model (Yule process).

With a pure birth model, the rate of only one event (birth of two new nodes from a parent one) is required. This model can be used for the simulation of species trees where the birth rate can be seen as the speciation rate. From a time \\( t_0 )\\ corresponding to the starting time of evolution with a single ancestral species, new branches leading to new species are generated. This process is also referred to as the Yule process [1]. 

This is a really simple and quite unrealistic model since it entirely ignore the dynamics of evolution. At any time \\( t_i \\) there are \\( n_i \\) number of individuals in the population, and each individual is able to independently give birth to an offspring at constant rate \\( \sigma \\). At most one birth can occur within a time interval  \\( ]t_i, t_i + \Delta t[ \\) and after a birth has happened, the parent and child evolve independently.

<div id="yule-process">
</div>

For a straightforward implementation (code below), we start with a pool of extant species of size one. After a waiting time which depends upon the speciation rate \\( \sigma \\), this single ancestral species is replaced by its two children in the pool of extant species. The process is then repeated several times, with a random extant species, until a stopping criterion is reached. The most common stopping criteria are : 

- Total evolution time
- Total number of extant species (number of leaves)

See <a href="#file-tree_simul-py-L19">birth_only_tree</a> at the end of the post for an implementation in python that use <strong><a href="http://etetoolkit.org">ETE</a></strong> and can be easily adapted to your need.

And an example of output is :

<pre>(((T1:0.523488,T2:0.523488)1:0.422493,(T3:0.79018,(T4:0.547927,T5:0.547927)1:0.242252)1:0.155801)1:0.294638,(T6:1.12469,T7:1.12469)1:0.115934);</pre>

![Pure Yule tree]({{ site.baseurl }}/public/images/tree_simul/yule.svg "Ultrametric tree")


As you can see, the simulated tree is perfectly ultrametric.

## Simulating a tree in a birth-death model

In a birth-death model, lineage can go extinct. Therefore, a rate for death events should additionally be be considered. Similar to the Yule process, in a small interval of time \\( ]t_i, t_i + \Delta t[ \\) , only one of the three following transitions are possible: `one death`, `one birth` or `no event`. 
The birth-death model is one of the most used for the simulation of species trees. Implementation of this model (<a href="#file-tree_simul-py-L84">birth_death_tree</a> ) follows a structure similar to the one above. The most important difference is to ensure that extinct species cannot be selected to give birth and should not have their branch lengths updated either. 

And an example of output is :

<pre>(((T1:1.59667,T2:1.59667)1:0.63352,(T3:1.66694,T4:1.66694)1:0.563248)1:1.12677,(:1.27368,:0.239961)1:0.268062);</pre>

![Birth-death tree]({{ site.baseurl }}/public/images/tree_simul/bd.svg "Birth-death tree")

The lost nodes are shown in gray dashed line. You can see that the corresponding lineage are at a shorter distance from the root compared to extant species. 

In the second part of this post, I will show how to simulate a gene family evolution  with duplication, loss and transfer events, given a species tree. 

## Implementation


{% gist 09f2b2f47713da9d465a1a7f72078fe7 %}


### References

[1] Yule G., “<em>A mathematical theory of evolution, based on the conclusions of Dr. J.C. Willis, F.R.S.</em>”, Philosophical Transactions of the Royal Society of London, Series B, vol. 213, pp. 21–84, 1924

<style type="text/css">

.node circle {
  fill: #999;
}

.link {
  fill: none;
  stroke: #555;
  stroke-opacity: 0.4;
  stroke-width: 1.5px;
}
</style>

<script src="http://d3js.org/d3.v4.min.js"></script>
<script type="text/javascript">

// set the dimensions and margins of the diagram
var margin = {top: 40, right: 50, bottom: 40, left: 50},
    width = 550 - margin.left - margin.right,
    height = 350 - margin.top - margin.bottom;

var sfid = 10, i = 0;
var root = {id:i, name:"root", parent:null},
    leaves = [root],
    tree = d3.cluster().size([width, height]);
    duration = 1200,
    timer = setInterval(update, duration);

var svg = d3.select("#yule-process").append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom),
    g = svg.append("g")
      .attr("transform",
            "translate(" + margin.left + "," + margin.top + ")");

var nodes = d3.hierarchy(root);
    nodes = tree(nodes);

root.x0 = width /2;
root.y0 = 0

g.selectAll("circle")
    .data(tree(nodes))
  .enter().append("svg:circle")
    .attr("class", "node")
    .attr("r", 10)
    .attr("cx", x)
    .attr("cy", y);

function update() {
  if (leaves.length >= 5){
    root = {id:i, name:"root", parent:null};
    leaves = [root];
    data = [root];
    nodes = d3.hierarchy(root);
    nodes = tree(nodes);
root.x0 = width /2;
root.y0 = 0
  g.selectAll("circle").remove()

  g.selectAll("path.link").remove();

g.selectAll("circle")
    .data(tree(nodes))
  .enter().append("svg:circle")
    .attr("class", "node")
    .attr("r", 10)
    .attr("cx", x)
    .attr("cy", y);

  }
  // Add a new datum to a random parent.
  // Add a new datum to a random parent.
  var source_pos = ~~(Math.random() * leaves.length),   source = leaves[source_pos];
  
  var d1 = {id:++i, parent:source} , d2 ={parent:source, id:++i};
  
  source.children = [d1, d2]
  leaves.splice(source_pos,1);
  leaves.push(d1);
  leaves.push(d2);

  
  // Compute the new tree layout. We'll stash the old layout in the data.
  nodes = d3.hierarchy(root);
  tree(nodes)
var node = g.selectAll(".node")
    .data(nodes.descendants(), nodeId)
  
node.enter().append("svg:circle")
      .attr("class", "node")

      .attr("r", 10)
      .attr("cx", function(d) { if(d.parent) return d.parent.data.x0; else return root.x0; })
      .attr("cy", function(d) { if(d.parent) return d.parent.data.y0; else return root.y0; })
  .attr("fill", "blue")
    .transition()

      .duration(1000)
      .attr("fill", "black")
      .attr("cx", x)
      .attr("cy", y);

  node.exit().remove();

  // Transition nodes to their new position.
  node.transition()
      .duration(1000)
      .attr("cx", x)
      .attr("cy", y);

// adds each node as a group
var link = g.selectAll("path.link")
      .data(nodes.links());

  // Enter any new links at the parent's previous position.
  link.enter().insert("path", "circle")
      .attr("class", "link")
      .attr("d", function(d) {
        var o = {x: d.source.data.x0, y: d.source.data.y0};
        return diagonal(o, o);      
      })
    .transition()
      .duration(1000)
      .attr("d", function(d) { return diagonal(d.source, d.target);})

  // Transition links to their new position.
  link.transition()
      .duration(1000)
      .attr("d", function(d) { return diagonal(d.source, d.target);})

}


function diagonal(d, p) {

  return "M" + d.x + "," + d.y
         + "C" + d.x + "," + (d.y + p.y) / 2
         + " " + p.x + "," +  (d.y + p.y) / 2
         + " " + p.x + "," + p.y;
  }

function x(d) {
  return d.data.x0 = d.x;
}

function y(d) {
  return d.data.y0 = d.y;
}

function nodeId(d) {
  return d.data.id;
}


</script>
