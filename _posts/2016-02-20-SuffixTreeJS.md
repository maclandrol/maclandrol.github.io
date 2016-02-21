---
layout: post
title:  "SuffixTreeJS"
date:   2015-12-22
tags: [en, js]
---

Ukkonen's Algorithm for generalized suffix tree in javascript, using d3js for visualization : [SuffixTreeJS](http://mrnoutahi.com/SuffixTreeJS/)
<!--more-->
This post will be short. Last year, I was a teaching assistant for a course about "Algorithms in Bioinformatics". Some of the students ask me for a javascript code to build generalized suffix trees in order to practice. After a not-so-extensive search, I only found this : [http://visualgo.net/suffixtree.html](http://visualgo.net/suffixtree.html), which do not build generalized suffix tree. 

As lazy as I was, I just made this : https://jsfiddle.net/tgqL2444/3/ using [Word Trees from Google chart](https://developers.google.com/chart/interactive/docs/gallery/wordtree) and send it to them. Fast forward a couple of days ago, I found this little snippet on d3js website:  [Interactive d3.js tree diagram](http://bl.ocks.org/d3noob/8375092) and thought to myself, maybe I can just use d3js to build a small interactive _generalized suffix tree__ viewer, and here it is: [SuffixTreeJS](http://mrnoutahi.com/SuffixTreeJS/).  It does the job but can still be improved.
