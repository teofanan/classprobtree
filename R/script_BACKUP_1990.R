# Copyright Jelena Vasic
#
#' Calculate the entropy of a vector
#'
#' \code{entropy} uses Claude Shannon's entropy formula
#'
#' @param x the vector
#' @return The entropy value
#' @export
#' @seealso \code{\link{information_gain}}
entropy<-function(x){
   if (class(x)!="factor") x<-factor(x)
   e<-0
   if (length(x)!=0) {
      for (curr_level in levels(x)) {
         sub_len<-length(x[x==curr_level])
         if (sub_len!=0) {
            p<-length(x[x==curr_level])/length(x)
            e<-e-p*log(p, 2)
         }
      }
   }
   e
}

#' Calculate the Information Gain for a set split
#'
#' \code{information_gain} calculates the Information Gain resulting from the splitting of a set. This is the difference between the entropy of the entire set and the weighted entropies of the
#'
#' @param x the target vector i.e. the entire set of labels (target variable values)
#' @param y the splitting vector i.e. the vector of values for the splitting attribute
#' @return The calculated Information Gain (IG)
#' @export
#' @seealso \code{\link{entropy}} \code{\link{best_information_gain}}
information_gain<-function(x, y) {
   if (class(x)!="factor") x<-factor(x)
   if (class(y)!="factor") y<-factor(y)
   if (length(x)!=length(y)) {
      paste("ERROR: Information Gain cannot be calculated because attribute and target vector lengths differ.")
      return
   }
   ig<-entropy(x)
   if (length(x)!=0) {
      for (curr_att_level in levels(y)) {
         sub_x<-x[y==curr_att_level]
         ig<-ig-entropy(sub_x)*length(sub_x)/length(x)
      }
   }
   ig
}

#' Find the attribute that provides the best information Information Gain as a splitting attribute
#'
#' \code{best_information_gain} investigates a set of attributes of a dataset to find the one that provides the best Information Gain as a splitter for a given target variable.
#'
#' @param x the target vector i.e. the entire set of labels (target variable values)
#' @param f a data frame with the candidate attributes as columns
#' @return A list with two members: the value of the best information gain and the column index in \code{f} of the 'winning' attribute
#' @export
#' @seealso \code{\link{information_gain}}
best_information_gain<-function(x, f) {
   best_ig<-0
   best_attrib<-0
   i<-0
   for (col in f) {
      i<-i+1
      curr_ig<-information_gain(x, col)
      if (best_attrib==0 | curr_ig>best_ig) {
         best_ig<-curr_ig
         best_attrib<-i
      }
   }
   list(best_ig, best_attrib)
}

#' Build a classification tree
#'
#' \code{classification_tree} builds a classification tree using Information Gain as the splitting criterion
#'
#' @param df the data frame containing the dataset
#' @param target the target variable, specified either using a character string or the column index in df
#' @param nxt the method of growing the tree, which should have one of the following two values: \code{"seq"} (sequentially, level by level) and \code{"ig"} (by highest IG, node by node)
#' @param stopc the value at which three building should stop, representing the \strong{number of levels} when \code{nxt="seq"} and the \strong{number of nodes} when \code{nxt="ig"}
#' @return A dataframe containing a row for each tree node, with the following columns:
#'
#'       \strong{Parent}: row index of parent node
#'
#'       \strong{Value}:  the value of the parent's splitting attribute that defines this node
#'
#'       \strong{Attrib}: the splitting attribute for this node
#'
#'       \strong{Level}:  the level of the node (starts with 1 for the root node, has value 2 for the root nodes' child nodes etc.)
#'
#'       \strong{IG}:     the Information Gain achieved by the splitting attribute at this node
#'
#'       \strong{Split}:  a logical value indicating whether the node has been split or not (some nodes may have a splitting attribute assigned but have IG=0, resulting in them not being split)
#'
#'       \strong{P_<cls>}: proportion of the target class in the subset represented by the node, one for each of the target classes
#'
#' @export
#' @seealso \code{\link{plot_class_tree}}
classification_tree<-function(df, target, nxt="seq", stopc=1000) {

   # convert the target to an index in the df
   if (class(target)=="character") {
      if (length(which(names(df)==target))!=0) {
         target<-which(names(df)==target)[1]
      } else {
         cat("ERROR: The target variable name is invalid (", target, ")\n")
         return
      }
   } else {
      if (class(target)=="numeric") {
         target<-as.integer(target)
      } else if (class(target)!="integer") {
         cat("ERROR: The target variable specification is of an invalid type (", class(target), ")\n")
         return
      }
      if (target<1|target>length(df)) {
         cat("ERROR: The target variable index is out of bounds (", target, ")\n")
         return
      }
   }

   # prepare the tree, in the form of a dataframe that stores the nodes as rows; the first row is the root node
   tree<-data.frame(Parent=NA, Value=NA, Attrib=NA, Level=1, IG=NA, Split=FALSE)
   index<-1
   tv<-df[[target]]
   for(cls in levels(tv)) {
      tree<-cbind(tree, data.frame(length(tv[tv==cls])/length(tv)))
      names(tree)[length(names(tree))]<-cls
   }


   # function that gets the attributes unused by nodes above n in the hierarchy of tree t
   available_attributes<-function(n) {
      used_attributes<-target
      n<-tree[n,"Parent"]
      while (!is.na(n)) {
         used_attributes<-c(used_attributes, tree[n,"Attrib"])
         n<-tree[n,"Parent"]
      }
      c(1:length(df))[-used_attributes]
   }

   # function that returns a vector with the subset indices for a particular node n in tree t
   subset_selector<-function(n) {
      # start with an all-TRUE selector vector
      selector<-rep(TRUE, nrow(df))
      # then iterate through the parents of the current node, at each one selecting only the instances
      # that have the attribute values associated with the node (the value will be stored in the child
      # and the splitting attribute name in the parent)
      while (n>1) {
         selector<-selector & df[[tree[tree[n, "Parent"], "Attrib"]]]==tree[n, "Value"]
         n<-tree[n, "Parent"]
      }
      selector
   }

   # returns the subset using the subset selector
   # TODO cacheing for better performance, if needed
   get_subset<-function(n) {
      df[subset_selector(n),]
   }

   # returns indicator as to whether an attribute has been identified that produces the best information gain
   set_information_gain<-function(n) {
      df_subset<-get_subset(n)
      curr_available_attribs<-available_attributes(n)
      result<-best_information_gain(df_subset[[target]], df_subset[curr_available_attribs])
      tree[n, "IG"]<<-result[[1]]
      if (result[[2]]==0) {
         tree[n,"Attrib"]<<-0
      } else {
         tree[n,"Attrib"]<<-curr_available_attribs[result[[2]]]
      }
      tree[n, "Attrib"]!=0
   }


   # function that return the next node to be split in tree t based on best information gain value
   next_to_split_by_ig<-function() {
      # go through all the nodes that don't have IG calculated and calculate+set it for them
      for (n in which(is.na(tree[,"IG"]))) {
         set_information_gain(n)
      }
      best_ig<-0
      best_node<-0
      # find the best IG among unsplit nodes
      for (n in which(!tree[,"Split"])) {
         if(best_node==0 | tree[n, "IG"]>best_ig) {
            best_ig<-tree[n, "IG"]
            best_node<-n
         }
      }
      best_node
   }

   split_set<-function(n) {
      # continue only if an attribute has been found and the information gain is positive
      if (tree[n,"IG"]!=0) {
         # set the parent as split
         tree[n, "Split"]<<-TRUE
         for (val in levels(df[, tree[n, "Attrib"]])) {
            newrow<-data.frame(Parent=n, Value=val, Attrib=NA, Level=tree[n, "Level"]+1, IG=NA, Split=FALSE)
            # add placeholder values for the probabilities
            classes<-levels(df[[target]])
            for(cls in classes) { newrow<-cbind(newrow, 0.0) }
            names(newrow)<-names(tree)
            tree<<-rbind(tree, newrow)

            # add the probabilities
            row_index<-nrow(tree)
            tv_subset<-get_subset(row_index)[[target]]
            column_index<-length(tree) - length(classes) + 1
            for(cls in classes) {
               tree[row_index, column_index]<<-length(tv_subset[tv_subset==cls])/length(tv_subset)
               column_index<-column_index+1
            }
         }
         TRUE
      } else {
         FALSE
      }
   }

   if (nxt=="seq") {
      n<-1
      # iterate, looking for best splitter, until all attributes used up
      while(n<=nrow(tree) & tree[n, "Level"]<stopc) {
         if (set_information_gain(n)) {
            split_set(n)
         }
         n<-n+1
      }
   }
   else if (nxt=="ig") {
      count<-0
      while(count<stopc) {
         n<-next_to_split_by_ig()
         # break if all the nodes have been split
         if (n==0) {
            break
         }
         # break if the set hasn't been split based on the last found n - this means that the IG is 0
         if (!split_set(n)) {
            break
         }
         count<-count+1
      }
   }
   tree
}

#' Plots a classification tree
#'
#' \code{plot_class_tree} plots a classification tree from data in the format created by \code{\link{classification_tree}}.
#'
#' @param tree the tree in the form of a dataframe like that one created by function \code{\link{classification_tree}}
#' @param attrib_names the names of the attributes from the dataset (the tree dataframe stores these as indices and the names make for a more readable plot)
#' @param fill_page fill the page with the diagram
#'
#' @export
#' @seealso \code{\link{classification_tree}}
plot_class_tree<-function(tree, attrib_names, fill_page=FALSE, probabilities=FALSE) {
   hpositions<-rep(-1, nrow(tree))
   # this is a recursive function; as it is anonymous it is referenced using sys.function(0)
   hpositions<-(function(hpositions, node){
      child_nodes<-which(tree[,"Parent"]==node)
      if (length(child_nodes) > 0) {
         for (child in child_nodes) {
            hpositions<-sys.function(0)(hpositions, child)
         }
         hpositions[node]<-(hpositions[child_nodes[1]]+hpositions[child_nodes[length(child_nodes)]])/2
      } else {
         hpositions[node]<-max(hpositions)+1
      }
      hpositions
   })(hpositions, 1)

   depth<-(function(node){
      curr_depth<-tree[node,"Level"]
      for (child in which(tree[,"Parent"]==node)) {
         child_depth<-sys.function(0)(child)
         if (child_depth>curr_depth) {
            curr_depth<-child_depth
         }
      }
      curr_depth
   })(1)


   # configuration
   vunit_multiplier<-5
   hunit_multiplier<-1
   text_halfspacing<-0.7

   # find the maximal position
   par(mar=c(0,0,0,0))
   if(probabilities) {
      # rough adjustment - TODO could be made more precise
      plot(c(-0.2,1.2), c(-0.3,1.1), type='n')
   } else {
      plot(c(-0.1,1.1), c(-0.1,1.1), type='n')
   }

   text_height<-strheight("M")
   fitting_text_height<-1/(depth-1)/vunit_multiplier

   if (fill_page | text_height>fitting_text_height) {
      par(cex=fitting_text_height/text_height)
      text_height<-fitting_text_height
   }

   vunit<-text_height*vunit_multiplier

   if (fill_page | max(hpositions)*hunit_multiplier*vunit>1) {
      hunit_multiplier<-1/vunit/max(hpositions)
   }

   hunit<-vunit*hunit_multiplier
   vertical_offset_from_char_centre<-text_height*text_halfspacing
   divline_halflen<-text_height*1.5
   offset_x<-0.5-max(hpositions)*hunit/2
   offset_y<-0.5-(depth-1)*vunit/2
   poz_x<-function(x) { x*hunit+offset_x }
   poz_y<-function(y) { y*vunit+offset_y }

   for (node in 1:nrow(tree)) {
      # value
      node_x<-poz_x(hpositions[node])
      toptext_y<-poz_y(depth-tree[node,"Level"])
      if (node==1) {
         text(label="root", x=node_x, y=toptext_y)
      } else {
         text(label=tree[node,"Value"], x=node_x, y=toptext_y)
      }

      # child nodes
      child_nodes<-which(tree[,"Parent"]==node)
      if (length(child_nodes) > 0) {
         # underline
         divline_y<-toptext_y-vertical_offset_from_char_centre
         lines(c(node_x-divline_halflen, node_x+divline_halflen), c(divline_y, divline_y), lwd=2)

         # splitting attribute
         bottomtext_y<-divline_y-vertical_offset_from_char_centre
         text(label=attrib_names[tree[node,"Attrib"]], x=node_x, y=bottomtext_y)

         # tree branches
         horizontal_branchline_y<-bottomtext_y-vertical_offset_from_char_centre
         horizontal_branchline_left_x<-poz_x(hpositions[child_nodes[1]])
         horizontal_branchline_right_x<-poz_x(hpositions[child_nodes[length(child_nodes)]])
         lines(c(horizontal_branchline_left_x, horizontal_branchline_right_x), c(horizontal_branchline_y, horizontal_branchline_y))
         for (child in child_nodes) {
            vertical_branchline_bottom_y<-toptext_y-vunit+vertical_offset_from_char_centre
            vertical_branchline_x<-poz_x(hpositions[child])
            lines(c(vertical_branchline_x, vertical_branchline_x), c(horizontal_branchline_y, vertical_branchline_bottom_y))
         }
      } else {
         # probability table
         if (probabilities) {
            # underline
            divline_y<-toptext_y-vertical_offset_from_char_centre
            lines(c(node_x-divline_halflen, node_x+divline_halflen), c(divline_y, divline_y), lwd=2)

            probability_y<-divline_y-vertical_offset_from_char_centre
            for (cls in names(tree)[7:length(names(tree))]) {
               text(label=paste(cls, ": ", tree[node,cls], sep=""), x=node_x, y=probability_y)
               probability_y<-probability_y-vertical_offset_from_char_centre*2
            }
         }
      }
   }

   par(cex=1)
}


#' Get probabilities of class membership for a set of data instances
#'
#' \code{class_membership_probabilities} returns the probabilities of a set of given instances belonging to a specified target class.
#'
#' @param tree the tree with classification and probability estimation information, built using the function \code{\link{classification_tree}}
#' @param df the dataframe containing the data instances
#' @param target_class the target variable's value to be used as the class for which membership probability is sought
#' @return A vector of probabilities.
#'
#' @export
#' @seealso \code{\link{rank_instances}}
class_membership_probabilities<-function(tree, df, target_class) {

   probability_vector<-c()
   for (i in 1:nrow(df)) {
      n<-1
      while (tree[n,"Split"]) {
         n<-which(tree[["Value"]]==df[i,tree[n,"Attrib"]] & tree[["Parent"]]==n)
      }
      probability_vector<-c(probability_vector, tree[n, target_class])
   }
   probability_vector
}


#' Ranks instances in a data frame using a classification tree
#'
#' \code{rank_instances} ranks the instances in a data set by probability of them belonging to a specified target class.
#'
#' @param tree the tree with classification and probability estimation information, built using the function \code{\link{classification_tree}}
#' @param df the dataframe containing the data instances to be ranked
#' @param target_class the target variable's value to be used as the class for ordering (i.e. the higher the probability of the instance being in that class, the higher its rank in the list)
#' @return A vector of instance indices in df, ordered by decreasing probability of membership of class \code{top_class}.
#'
#' @export
#' @seealso \code{\link{classification_tree}}
rank_instances<-function(tree, df, target_class) {
   c(1:nrow(df))[order(class_membership_probabilities(tree, df, target_class), decreasing=TRUE)]
}
