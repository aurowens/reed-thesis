---
title: "Introduction to R and RStudio"
output:
  html_document:
    highlight: pygments
    theme: cerulean
    toc: true
    toc_float: true
---

# Literature Review

* * *

##*Statistical Modeling: The Two Cultures*

- Leo Breiman, 2001, *Statistical Science*,  Vol. 16, No. 3, 199–231

Breiman presents an argument that for statistics to grow as a field, we must move away from modeling under the assumption that the data arrives from a stochastic data model. 

In many ways this paper is *the paper* for this endeveour. Breiman identifies the singular issue with inferential statistics using tree-based methods: variable importance and inferential techniques from logistic regression seldom agree. 

* * *

##*Hierarchical Testing of Variable Importance*

- Nicolai Meinshausen

Meinhausen looks at "variable importance" in the context of highly correlated predictors in a linear model. He proposes clustering the variables and testing each with a p value adjusted for their place in the cluster. I have a few things that I am unsure on with this paper:

- Why not PCA? 
- Meinhausen cites the Brieman paper and uses the term "variable importance" but is almost certainly only refering to traditional significance tests.

While not directly germane to random forests, this paper presents an interesting question: maybe the miscommunication between variable importance and significance is partially caused by correlation in the variables?

* * *

##*Conditional Variable Importance for Random Forests*

- Carolin Strobl, Anne-Laure Boulesteix, Thomas Kneib, Thomas Augustin and Achim Zeileis, Technical Report Number 23, 2008, University of Munich

Strobl et al respond to Breiman, Cutler's permutation method for identifying variable importance by providing a new method that tests for conditional independence. 

I have a few questions about this paper:

- Ultimately, their findings rest on mtry being a "sufficiently large" value. There still seems to be some arbitrariness in that specification.

- Does this method produce variable importance plots that "agree" with other models? 

This paper seems terribly releavent to my thesis as it addresses and makes clear some of the things that Breiman mentions in the 2001 paper. 

* * *

##*Random forests-classification description*

-Breiman and Cutler 2007

This is how Strobl cites the method Breiman and Cutler developed for a permuted variable importance significance test. I haven't been able to find this as a paper, and Strobl's bibliography suggests this was a website in 2007. Breiman and Cutler's current website https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm#micro2, includes a small paragraph about their method.  

* * *

##*Variable importance in binary regression trees and forests*

- Hemant Ishwaran, Electronic Journal of Statistics, Vol. 1 (2007) 519–537

Ishwaran develops another method for variable importance based on the assuption that predictors that split closest to the original node are more "important". Ishwaran's method follows from the intution of variable importance on a single tree to a theoritically - backed method for a RF. 

This paper is useful for the same reason as Sandri and Zuccolotto, as it provides another method for variable importance. 


* * *

##*Variable Selection Using Random Forests*

- Marco Sandri and Paola Zuccolotto, Dipartimento Metodi Quantitativi, Universit`a di Brescia

Sandri and Zuccolotto note four different ways to measure variable importance from a RF. The first one is Breiman and Cutler's permutation stradegy. They then propose a new method of variable selection that combines all four:

1. Create a centroid with coordinates given by the average (or median? *they included this note but I am not sure I follow it*) of the four measures 

2. Calculate the distance of each point variable from the centroid and arrange in descending order

3. Choose a threshold m, and say that all variables with distance > m are "important"

There is an analougous method involving PCA the authors briefly mention as well, since the values of variable importance recieved from the four methods are most often correlated. 

This is relevent to my thesis as it is potentially another method of variable importance to analyze.

* * *

##Works Cited by these Authors to Explore:

- *Classification with Random Forests: the theoretical framework*,  Sandri and Zuccolotto (2004)

- 
