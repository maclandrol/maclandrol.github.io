---
layout: post
title:  "Fisher's Exact test for mxn contingency table"
date:   2015-12-31
tag: en, stat
---

I was working on a project where I needed to compare categorical data in other to determine if there is any association between them. This could be done with a chi2 test for independence which has implementation in most programming languages. 
<!--more-->
The problem was that I had low count in my contingency table. Although it's commonly accepted that the Cochran's rule **of at least 5 count in each cell of the expected table is too strict**, I actually had sometimes 0 in my observed table and really low total count. Thus, using an exact test is probably more suited. In fact, [McDonald recommend, as a rule of thumb, to use an exact test when the total count is lower than 1000](http://www.biostathandbook.com/small.html). 

So I decided to perform a Fisher's exact test, despite its controversies (see [Frank Harrell's post with detailled references on stats.stackexchange.com](http://stats.stackexchange.com/questions/14226/given-the-power-of-computers-these-days-is-there-ever-a-reason-to-do-a-chi-squa/14230#14230)). Unfortunately, I didn't find any implementation of the Fisher's exact test for larger than 2x2 table in python. To my knowledge, only [R offer a Fisher's exact test for table larger than 2x2](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/fisher.test.html). I then decided to rewrite the R function in python. How difficult would that been right ?

It turns out that the R function actually use a heavily edited C version of the FORTRAN subroutine FEXACT ([algorithm 643 by Mehta and Patel](http://dl.acm.org/citation.cfm?id=214326&picked=formats&preflayout=tabs)). My plan was to first find the [fortran 90 source code of FEXACT](http://jblevins.org/mirror/amiller/) and then use [f2py](http://docs.scipy.org/doc/numpy-dev/f2py/) to import the fortran module in my python code. I quickly learned that although this sound simple, it was nearly impossible without any basic notion in fortran. 

Since f2py couldn't handle real precision when ```SELECTED_REAL_KIND``` is used, I first created a [.f2py_f2cmap](https://sysbio.ioc.ee/projects/f2py2e/FAQ.html#q-what-if-fortran-90-code-uses-type-spec-kind-kind) file with ```{'real':{'dp':'double'}}``` in it in order to map any instance of dp to double. This wasn't enough, so I changed the structure of the FEXACT module :

```
MODULE Types

IMPLICIT NONE
INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(14, 60)

END MODULE

MODULE Fisher_Exact
USE Types
CONTAINS
...
END MODULE Fisher_Exact