---
layout: post
title:  "SuffixTreeJS"
date:   2016-02-20
tags: [en, js]
comments: true
categories: [Algorithm, Strings]
---

Ukkonen's Algorithm for generalized suffix tree in javascript, using d3js for visualization : [SuffixTreeJS](http://www.mrnoutahi.com/SuffixTreeJS/)

This post will be short. 

Last year, I was a teaching assistant for a course about "Algorithms in Bioinformatics". Some of the students ask me for a javascript code to build generalized suffix trees in order to practice. After a not-so-extensive search, I only found this : [http://visualgo.net/suffixtree.html](http://visualgo.net/suffixtree.html), which does not build a generalized suffix tree. 

As lazy as I was, I just made [this](https://jsfiddle.net/kgbhcpn7/) using [Word Trees from Google chart](https://developers.google.com/chart/interactive/docs/gallery/wordtree) and send it to them. Fast forward a couple of days ago, I found this snippet on d3js website:  [Interactive d3.js tree diagram](http://bl.ocks.org/d3noob/8375092) and thought it could be used to build a small interactive **_generalized suffix tree_** viewer, and here it is: [SuffixTreeJS](http://www.mrnoutahi.com/SuffixTreeJS/).  

It does the job but it can still be improved, so suggestions are welcome.
