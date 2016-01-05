---
layout: post
title:  "Fisher's Exact test for mxn contingency table"
date:   2016-01-03
tag: stat, python
---


I was working on a project where I needed to compare categorical data in other to determine if there is any association between them. For my particular problem, a chi2 test wouldn't work, so I needed a Fisher's exact test. Since there aren't any code in python to perform the Fisher's exact test for larger than 2x2 table, I decided to write my own. You can find it here : [FisherExact](https://github.com/maclandrol/FisherExact).
<!--more-->

Usually, you would use a chi2 test for independence (which is implemented in most programming languages) to compare association between categorical data. The problem was that I had low count in my contingency table.

Although it's commonly accepted that the Cochran's rule **of at least 5 count in each cell of the expected table is too strict**, I actually had really low total counts and cells with 0 count in my observed table. Thus, using an exact test is probably more suited. In fact, [McDonald recommend, as a rule of thumb, to use an exact test when the total count is lower than 1000](http://www.biostathandbook.com/small.html). 

So I decided to perform a Fisher's exact test, despite its controversies (see [Frank Harrell's post with detailled references on stats.stackexchange.com](http://stats.stackexchange.com/questions/14226/given-the-power-of-computers-these-days-is-there-ever-a-reason-to-do-a-chi-squa/14230#14230)). Unfortunately, I didn't find any implementation of the Fisher's exact test for larger than 2x2 table in python. To my knowledge, only [R offer a Fisher's exact test for table larger than 2x2](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/fisher.test.html). Since the main part of my project was already written in python, I decided to rewrite the R function in python for consistency. How difficult would have that been ? right ?

It turns out that the R function actually use a heavily edited C version of the FORTRAN subroutine FEXACT ([algorithm 643 by Mehta and Patel](http://dl.acm.org/citation.cfm?id=214326&picked=formats&preflayout=tabs)). My plan was to first find the [fortran 90 source code of FEXACT](http://jblevins.org/mirror/amiller/) and then use [f2py](http://docs.scipy.org/doc/numpy-dev/f2py/) to import the fortran module in my python code. I quickly learned that, although this sound simple, it was nearly impossible without any basic notion in fortran. 

Since f2py couldn't handle real precision when ```SELECTED_REAL_KIND``` is used, I first created a [.f2py_f2cmap](https://sysbio.ioc.ee/projects/f2py2e/FAQ.html#q-what-if-fortran-90-code-uses-type-spec-kind-kind) file with ```{'real':{'dp':'double'}}``` in it in order to map any instance of ```dp``` to double. This wasn't enough, so I changed the structure of the FEXACT module to this :

{% highlight fortran %}
MODULE Types

IMPLICIT NONE
INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(14, 60)

END MODULE

MODULE Fisher_Exact
USE Types
CONTAINS
...
END MODULE Fisher_Exact

{% endhighlight %}

I also made the input and output more explicit, using ```INTENT(OUT)``` and ```INTENT(IN)```.

I wrote the python code based on the R version, so it should essentially do the same thing. For 2x2 contingency  table, the [fisher_exact](http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.stats.fisher_exact.html) function from scipy.stats is used. Since R also provide p-values computation by Monte Carlo simulation, I tried to offer the same thing in my knockoff version. A quick glance at R sources, reveal that the fisher_exact simulation use a C version of the [rcont2 subroutine](http://people.sc.fsu.edu/~jburkardt/f_src/asa159/asa159.html) (written in Fortran90). After some struggle, during which, I learn about [array indexing order in C vs Fortran](http://docs.scipy.org/doc/numpy-1.10.0/reference/internals.html#multidimensional-array-indexing-order-issues), in-place modification of argument with ```INTENT(INPLACE)``` and many more, I was able to make it work. 

I also added mid-P correction (see [this document](http://www.biddle.com/documents/AIToolkit_Important_08-2010.pdf), [this one](http://www.sheffield.ac.uk/polopoly_fs/1.43998!/file/tutorial-9-fishers.pdf)and [Lydersen, Fagerland and Laake, 2009](http://onlinelibrary.wiley.com/doi/10.1002/sim.3531/epdf) for reference), which is not present in R. 

The source code is available here : https://github.com/maclandrol/FisherExact.


### Comparing R's fisher.test to my fisher_exact.


<table>
  <thead>
    <tr>
      <th>Data</th>
      <th>Parameters</th>
      <th>pval : fisher.test (R)</th>
      <th>pval : FisherExact (python)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>
	<table border="0">         
           <td>3</td>
           <td>1</td>
           <tr>
           <td>1</td>
           <td>3</td>
           </tr>
           </table>
      </td>
      <td>default (exact)</td>
      <td>0.4857</td>
      <td>0.4857</td>
    </tr>
    <tr>
      <td>
	<table border="0">
           <td>1</td>
           <td>3</td>
           <td>10</td>
           <td>6</td>
           <tr>
           <td>2</td>
           <td>3</td>
           <td>10</td>
           <td>7</td>
           </tr>
			<tr>
           <td>1</td>
           <td>6</td>
           <td>14</td>
           <td>12</td>
           </tr>
			<tr>
           <td>0</td>
           <td>1</td>
           <td>9</td>
           <td>11</td>
           </tr>
           </table>
      </td>
      <td>
      	default (exact)
      </td>
      <td>
      0.7827
      </td>
      <td>
      0.78268
      </td>
    </tr>
     <tr>
      <td>
	<table border="0">
           <td>1</td>
           <td>3</td>
           <td>10</td>
           <td>6</td>
           <tr>
           <td>2</td>
           <td>3</td>
           <td>10</td>
           <td>7</td>
           </tr>
			<tr>
           <td>1</td>
           <td>6</td>
           <td>14</td>
           <td>12</td>
           </tr>
			<tr>
           <td>0</td>
           <td>1</td>
           <td>9</td>
           <td>11</td>
           </tr>
           </table>
      </td>
      <td>
      	simulated p-values, replicates=1e5
      </td>
      <td>
      0.7829
      </td>
      <td>
      0.7823
      </td>
    </tr>
         <tr>
      <td>
	<table border="0" CELLPADDING="1" CELLSPACING="1">
           <td>1</td>
           <td>3</td>
           <td>10</td>
           <td>6</td>
           <tr>
           <td>2</td>
           <td>3</td>
           <td>10</td>
           <td>7</td>
           </tr>
			<tr>
           <td>1</td>
           <td>6</td>
           <td>14</td>
           <td>12</td>
           </tr>
			<tr>
           <td>0</td>
           <td>1</td>
           <td>9</td>
           <td>11</td>
           </tr>
           </table>
      </td>
      <td>
      	hybrid
      </td>
      <td>
      0.7827
      </td>
      <td>
      0.78268
      </td>
    </tr>
  </tbody>
  </table>


As expected, we have the same p-values for 2x2 table (here the scipy fisher_exact is used) and for tables larger than 2x2 with the default parameters. For the simulated p-value, R doesn't offer a seed option, so it's difficult the two p-values, but I think it's working because, even after several runs, the results were always consistent.
In hybrid mode, an approximation based upon asymptotic chi-squared probabilities is used instead of Fisher exact test probabilities. You can set some arguments in FEXACT (```expect=5.0, percnt=80.0 and emin=1.0```), to obtain the **'Cochran'** condition. In R's source code, ```percnt=180.0``` is used instead : 
  
{% highlight R %}
else if(hybrid) {
	## Cochran condition for asym.chisq. decision:
	PVAL <- .Call(C_Fexact, x, c(5, 180, 1), workspace, mult)
	## ~ Added by me for explanation
	## ~ This line call Fexact from the C source code, on matrix x
} else {
	##  expect < 0 : exact
	PVAL <- .Call(C_Fexact, x, c(-1, 100, 0), workspace, mult)
}
   
{% endhighlight %}

I'm not sure if this is an error, but the comment preceding said "Cochran condition", so I guess something is not right. [I report it as a bug](https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16654), although I'm not sure if this is a real bug or I'm just too stupid... Anyway, in my implementation I use ```percnt=80.0```, and unfortunately, on my test data, the hybrid and the exact mode return the same p-value.

In conclusion, you can now do Fisher's exact test on any mxn contingency table in python. I hope this will be useful to someone.