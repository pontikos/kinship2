
drop.columns <- function(d, columns) d[,-sapply(columns, function(i) which(colnames(d)==i))]

read.pedigree <- function(ped,kinship2=FALSE) {
    ped <- read.csv(ped, colClasses="character")
    rownames(ped) <- as.character(ped$ID)
    ped$Affection <- as.numeric(ped$Affection)
    if ('exclude' %in% colnames(ped)) ped <- ped[which(is.na(ped$exclude)),]
    ped2 <- pedigree(id=ped[,'ID'],dadid=ped[,'Father'],momid=ped[,'Mother'],aff=as.numeric(ped[,'Affection']),sex=as.numeric(ped[,'Gender']),missid="0")
    if (!'founder' %in% colnames(ped))  ped <- set.founder(ped)
    else ped$founder <- as.logical(ped$founder)
    if (kinship2) return(ped2kinship2(ped))
    else return(ped)
}

ped2kinship2 <- function(ped) {
    ped2 <- pedigree(id=ped[,'ID'],dadid=ped[,'Father'],momid=ped[,'Mother'],aff=as.numeric(ped[,'Affection']),sex=as.numeric(ped[,'Gender']),missid="0")
    return(ped2)
}

# takes a pedigree object as returned by read.pedigree
# and a list of variants with names corresponding to the individuals in the pedigree
plot.ped <- function(ped,variant=NULL,colors=NULL,legend.names=NULL,DNA=FALSE,AGE=FALSE,cex=.6,title='',output=NULL) {
    if (!is.null(variant)) {
        samples <- intersect(names(variant),rownames(ped))
        variant_ids <- c('variant')
        ped[samples,'variant'] <- variant[samples]
    } else {
        variant_ids <- c()
    }
    if (!is.null(colors)) {
        colors <- unlist(strsplit(colors, ','))
        if (length(colors)!=length(variant_ids)) stop('length of colors:',length(colors),'length of variant_ids:',length(variant_ids))
        names(colors) <- variant_ids
    }
    if (!is.null(legend.names)) {
        legend.names <- unlist(strsplit(legend.names, ','))
        if (length(legend.names)!=length(variant_ids)) stop('length of legend.names:',length(legend.names),'length of variant_ids:',length(variant_ids))
        names(legend.names) <- variant_ids
    }
    ped$box.col <- 'black'
    print(dim(ped2 <- ped2kinship2(ped)))
    print(length(ped$box.col))
    print(dim(ped2))
    cols <- c('0'='black','1'='darkgreen','2'='red')
    if (!is.null(output)) pdf(output,width=80)
    par(xpd=T)
    #a <- plot(ped2,cex=0.4,col=ped$box.col)
    a <- plot(ped2,cex=0.4)
    title(title)
    if (!is.null(variant)) {
        inc <- (cex/6)
        y <- a$y+inc
        for (variant_id in variant_ids) {
            y <- y+inc
            if (!is.null(colors)) {
                col <- colors[[variant_id]]
            } else {
                ped$col <- 'black'
                ped$col <- cols[as.character(ped[,variant_id])]
                ped$col[is.na(ped$col)] <- 'gray'
                col <- ped$col
            }
            text( a$x,y,label=ped[,variant_id],cex=cex,col=col)
        }
        #if (is.null(legend.names)) legend.names <- variant_ids
        if (!is.null(colors)) legend('topleft',legend.names, text.col=colors, bty='n')
        else legend('topleft',legend.names, bty='n')
        legend('bottomleft', c('wildtype=0','mutant=1'),bty='n')
    } else if (DNA) {
        inc <- cex/3
        y <- a$y+inc
        text( a$x,y,label=ped[,'DNA'],cex=opt$cex)
        legend('topleft','DNA',bty='n')
    } else if (AGE) {
        inc <- cex/3
        y <- a$y+inc
        text( a$x,y,label=ped[,'ImputedAge'],cex=opt$cex)
        legend('topleft','age',bty='n')
    }
    if (!is.null(output)) dev.off()
}


### GET
#the founders of the pedigree

get.parents <- function(ID, ped) {
    return(rbind(get.father(ID,ped),get.mother(ID,ped)))
}

get.father <- function(ID, ped) {
    return(ped[ped[ID,'Father'],])
}

get.mother <- function(ID, ped) {
    return(ped[ped[ID,'Mother'],])
}

get.spouse <- function(ID, ped) {
    if (is.na(ID)) return(NA)
    x <- setdiff(unique(ped[which(ped$Father==ID|ped$Mother==ID),c('Mother','Father')]),ID)
    if (nrow(x)==0) return(NULL)
    else return(ped[as.character(x),])
}

get.children <- function(ID, ped) {
    if (is.na(ID)) return(NA)
    children <- ped[which(ped$Father==ID | ped$Mother==ID),]
    if (nrow(children)>0) return(children) else return(NULL)
}

get.siblings <- function(ID, ped) {
    if (is.na(ID)) return(NA)
    rownames(ped) <- ped$ID
    ind <- ped[ID,]
    Father <- ind[['Father']]
    Mother <- ind[['Mother']]
    if (is.na(Father) | is.na(Mother)) return(NULL)
    siblings <- ped[which(ped$Father==Father & ped$Mother==Mother & ped$ID!=ID),]
    if (nrow(siblings)>0) return(siblings) else return(NULL)
}

get.youngest.sibling <- function(ID, ped) {
    if (is.na(ID)) return(NA)
    siblings <- get.siblings(ID,ped)
    return( siblings[which.min(siblings$ImputedAge),] )
}

get.oldest.sibling <- function(ID, ped) {
    if (is.na(ID)) return(NA)
    siblings <- get.siblings(ID,ped)
    return( siblings[which.max(siblings$ImputedAge),] )
}

get.subfamily <- function(ID, ped) {
    rownames(ped) <- ped$ID
    children <- get.children(ID,ped)
    if (nrow(children)==0) {
        return(ped[ID,])
    } else {
        if (!is.na(get.spouse(ID,ped))) ID <- c(ID,get.spouse(ID,ped))
        return(rbind(ped[ID,],do.call('rbind', lapply(children[,'ID'], function(child.id) get.subfamily(child.id, ped)))))
    }
}

get.generation <- function(ID, ped) {
    if (is.na(ID)) return(NA)
    rownames(ped) <- ped$ID
    ind <- ped[ID,]
    Father <- ind[['Father']]
    Mother <- ind[['Mother']]
    # if has no parents then this must be generation 0
    if (is.na(Father) & is.na(Mother)) return(0)
    else
    return(1+get.generation(Father,ped)+get.generation(Mother,ped))
}

get.cases <- function(ped,d=NULL) {
    if (is.null(d)) return(ped[which(ped$Affection==2),])
    else return(ped[which(ped$ID %in% intersect(ped[which(ped$Affection==2),'ID'],colnames(d))),])
}
    
get.controls <- function(ped,d=NULL) {
    if (is.null(d)) return(ped[which(ped$Affection==1),])
    else return(ped[which(ped$ID %in% intersect(colnames(d),ped[which(ped$Affection==1),'ID'])),])
}

get.founders <- function(ped, individuals=NULL) {
    if (!'founder' %in% colnames(ped)) ped <- set.founder(ped)
    if (is.null(individuals)) {
        return(ped[which(ped$founder),])
    } else {
        founders <- c()
        for (ind in individuals) {
            father <- ped[ind,'Father']
            mother <- ped[ind,'Mother']
            if (is.na(father)|is.na(mother)) founders <- unique(c(founders, ind))
            else founders <- unique(c(founders, get.founders(ped, father), get.founders(ped, mother)))
        }
    }
    return(unique(founders))
}

# returns the founder parent
get.founder <- function(ID, ped) {
    if (!'founder' %in% colnames(ped)) ped <- set.founder(ped)
    p <- get.parents(ID,ped)
    i <- which(p$founder)
    if (length(i)>0) return(p[i,])
    return(get.founder(p[1,'ID'],ped))
}

# returns the non-founder parent
get.non.founder <- function(ID, ped) {
    if (!'founder' %in% colnames(ped)) ped <- set.founder(ped)
    p <- get.parents(ID,ped)
    i <- which(!p$founder)
    if (length(i)>0) return(p[i,])
    return(get.non.founder(p[1,'ID'],ped))
}



#subfamily per founder: spouse+children
get.subfamilies <- function(ped) {
 subfamilies <- list()
 for (f in get.founders(ped)$ID) {
    spouse <- get.spouse(f,ped)
    if (!spouse %in% ped$ID) next
    subfam <- rbind(ped[f,], ped[spouse,])
    subfam <- rbind(subfam, get.children(f,ped))
    subfamilies[[f]] <- subfam
 }
 return(subfamilies)
}

# returns subfamily of individual by going down the tree
getSubFamily <- function(ID, ped) {
    x <- ped[which(ID==ped$ID),]
    if ( is.na(x['Subfamily']) ) {
        parent <- ped[which((x$Father==ped$ID | x$Mother==ped$ID) & !(ped$Father == 0 & ped$Mother ==0)),]
        # if no non-founder parent is found then look at spouse
        return(getSubFamily(parent$ID))
    } else {
        return(x['Subfamily'])
    }
}

# get all trios in pedigree
get.trios <- function(ped) {
    trios <- c()
    for (i in 1:nrow(ped)) {
        p <- ped[i,]
        if ( p['Father'] %in% ped[,'ID'] & p['Mother'] %in% ped[,'ID'] ) trios <- rbind(trios,p)
    }
    return(trios)
}

# get.root
get.root <- function(ped) {
  root <- data.frame()
  if (!'founder' %in% colnames(ped)) ped <- set.founder(ped)
  founders <- get.founders(ped)$ID
  for (f in founders) {
    spouse <- ped[f,'spouse']
    if (spouse %in% root$ID | f %in% root$ID) next
    if (as.logical(ped[spouse,'founder'])) root <- rbind(root,ped[c(f,spouse),])
  }
  return(root)
}

#
get.bits <- function(ped) {
    ped <- set.founder(ped)
    founders <- nrow(ped[ped$founder,])
    non.founders <- nrow(ped[!ped$founder,])
    return(2*non.founders - founders)
}

get.subfamilies <- function(ped, nbits=25)  {
    ped <- set.founder(ped)
    n <- nrow(ped)
    n.founders <- nrow(ped[ped$founder,])
    non.founders <- nrow(ped[!ped$founder,])
    max.non.founders <- (nbits+n.founders)/2
    # this might be too large
    choose(non.founders, max.non.founders)
}

### CHECKS

#check that both parents are in the pedigree
check.both.parents <- function(ped) {
    rownames(ped) <- ped$ID
    return( sapply( ped$ID, function(id) {
        ind <- ped[id,]
        Father <- ind[['Father']]
        Mother <- ind[['Mother']]
        if ( is.founder2(id,ped) ) return(TRUE)
        else return(  (Father %in% ped$ID) & (Mother %in% ped$ID) )
        }) )
}

is.founder2 <- function(ID,ped) {
    rownames(ped) <- ped$ID
    if (is.na(ID)) return(NA)
    ind <- ped[ID,]
    Father <- ind[['Father']]
    Mother <- ind[['Mother']]
    #return(is.na(Father) & is.na(Mother))
    return((!Father %in% ped$ID) & (!Mother %in% ped$ID))
}

check.is.founder <- function(ped) {
    rownames(ped) <- ped$ID
    return( sapply(ped$ID, function(id) is.founder2(id,ped)) )
}

#orphans have no parents AND are not parents to anyone
#return true or false
is.not.orphan <- function(id,ped) {
    ind <- ped[id,]
    return(!(is.na(ind$Father) & is.na(ind$Mother) & !(id %in% ped$Father |  id %in% ped$Mother)))
}
check.not.orphan <- function(ped) {
    rownames(ped) <- ped$ID
    return( sapply(ped$ID, function(id) is.not.orphan(id,ped)) )
}

#genders are consistent:
#males can only be fathers and females mothers
is.gender.consistent <- function(id,ped) {
    ind <- ped[id,]
    gender <- ind['Gender']
    return( !(gender==1 & (id %in% ind$Mother)) & !(gender==2 & (id %in% ind$Father)) )
}
check.gender.consistent <- function(ped) {
    rownames(ped) <- ped$ID
    return( sapply(ped$ID, function(id) is.gender.consistent(id,ped)) )
}

check.gender.known <- function(ped) {
    rownames(ped) <- ped$ID
    return( sapply(ped$ID, function(id) ped[id,'Gender'] %in% c(1,2)) )
}
    

# polygamy not supported
# every spouse must have a single spouse
# and the relation must be reciprocal
# NA means person is single
is.spouse.matched <- function(id,ped) {
    ind <- ped[id,]
    spouse <- get.spouse(id, ped)
    spouse.spouse <- get.spouse(spouse, ped)
    return(id==spouse.spouse)
}
check.spouse.matched <- function(ped) {
    rownames(ped) <- ped$ID
    return( sapply(ped$ID, function(id) is.spouse.matched(id,ped)) )
}

# all indviduals which have children should have a spouse
check.couples.children <- function(ped) {
    rownames(ped) <- ped$ID
    for (id in ped$ID) {
        spouse <- get.spouse(id,ped)
        n.children <- nrow(get.children(id,ped))
        if (n.children > 0 & is.na(spouse)) warning(id, 'children but no spouse')
        if (!is.na(spouse) & n.children==0) warning(id, 'spouse but no children')
    }
}
        
# kinship2 check
get.kinship2.pedigree <- function(ped) {
    return(pedigree(id=ped[,'ID'],dadid=ped[,'Father'],momid=ped[,'Mother'],aff=as.numeric(ped[,'Affection']),sex=as.numeric(ped[,'Gender']),missid="0"))
}

### SET

nrow2 <- function(x) {
    x <- nrow(x)
    if (is.null(x)) return(0)
    else return(x)
}

set.children <- function(ped) {
    ped$children <- sapply(ped$ID, function(id) nrow2(get.children(id,ped)))
    return(ped)
}

set.siblings <- function(ped) {
    ped$siblings <- sapply(ped$ID, function(id) nrow2(get.siblings(id,ped)))
    return(ped)
}

set.spouse <- function(ped) {
    ped$spouse <- sapply(ped$ID, function(id) get.spouse(id,ped))
    return(ped)
}

set.founder <- function(ped) {
    ped$founder <- sapply(ped$ID, function(id) is.founder2(id,ped))
    return(ped)
}

set.genereration <- function(ped) {
    rownames(ped) <- ped$ID
    ped$Generation <- NA
    ped$Generation <- sapply(ped$ID, function(id) get.generation(id,ped))
    ped$Generation <- sapply(ped$ID, function(id) if (is.founder2(id,ped)) return(ped[get.spouse(id,ped),'Generation']) else return(ped[id,'Generation']))
    return(ped)
}

# all descendants and spouses will be assigned to subfamily
setSubfamily <- function(ped,ID,Subfamily=NA) {
    ped[which(ped$ID==ID),'Subfamily'] <- Subfamily
    gender <- ped[which(ped$ID==ID),'Gender']
    if (gender==1) {
        (wife.id <- unique(ped[which(ped$Father==ID),'Mother']))
        ped[which(ped$ID==wife.id),'Subfamily'] <- Subfamily
    } else {
        (husband.id <- unique(ped[which(ped$Mother==ID),'Father']))
        ped[which(ped$ID==husband.id),'Subfamily'] <- Subfamily
    }
    children <- ped[which(ped$Father==ID|ped$Mother==ID),]
    for (i in seq_len(nrow(children))) {
        (child <- children[i,])
        ped[which(ped$ID==child['ID']),'Subfamily'] <- Subfamily
        ped <- setSubfamily(ped,child[['ID']],Subfamily)
    }
    return(ped)
}

#ped <- read.table('samples.ped',header=TRUE)
#
#ped$Subfamily <- NA
#ped <- setSubfamily(ped,3000,0)
#ped <- setSubfamily(ped,3001,1)
#ped <- setSubfamily(ped,3002,2)
#ped <- setSubfamily(ped,3003,3)
#ped <- setSubfamily(ped,3004,4)
#ped <- setSubfamily(ped,3005,5)
#ped <- setSubfamily(ped,3006,6) 
##Load pedigree
#seq.ped <- read.csv("DNA_pedigree_details.csv")
#seq.ped <- seq.ped[,c('ID','DNA','ImputedAge','Generation','uclex.sample')]
#print(dim(ped2 <- merge(ped, seq.ped, all.x=TRUE))) 
#write.csv(ped2,file='pedigree_details.csv',quote=FALSE,row.names=FALSE)

# mendelian error check
ME.check <- function(ped, d, fix=FALSE) {
    samples <- intersect(ped$ID,names(d))
    ped2 <- pedigree(id=ped[,'ID'],dadid=ped[,'Father'],momid=ped[,'Mother'],aff=as.numeric(ped[,'Affection']),sex=as.numeric(ped[,'Gender']),missid="0")
    errors <- c()
    for (i in 1:nrow(ped)) {
        id <- ped[i,'ID']
        father <- ped[i,'Father']
        mother <- ped[i,'Mother']
        if ( length(intersect(c(id,father,mother),samples))<3 ) next;
        i <- which(d[id] != 0  & d[father]==0 & d[mother]==0)
        if (length(i) > 0) {
            cat('ERR',sprintf("both parents %s %s are wt, child %s is not wt",mother,father,id), sep=':')
            cat('\n')
            errors <- c(errors, c(mother,father,id))
            if (fix) {
               #set child to wt
               #make trio NA
               d[id] <- NA
               d[father] <- NA
               d[mother] <- NA
            }
        }
        i <- which(d[id] == 0  & (d[father]==2 | d[mother]==2))
        if (length(i) > 0) {
            cat('ERR',sprintf("at least one parent %s %s is hom, child %s is wt",mother,father,id), sep=':')
            cat('\n')
            errors <- c(errors, c(mother,father,id))
            if (fix) {
                #make the hom parent het or make the child het?
               d[id] <- NA
               d[father] <- NA
               d[mother] <- NA
            }
        }
        i <- which(d[id] != 2  & d[father]==2 & d[mother]==2)
        if (length(i) > 0) {
            cat('ERR',sprintf("both parents %s %s are hom, child %s is not hom",mother,father,id), sep=':')
            cat('\n')
            errors <- c(errors, c(mother,father,id))
            if (fix) {
               #set child to hom
               d[id] <- NA
               d[father] <- NA
               d[mother] <- NA
            }
        }
    }
    print(errors <- unique(errors))
}

# not finished
random.walk <- function(ped,nbits=40,id=sample(ped$ID,1)) {
     if (is.founder2(ped=ped,ID=ID)) nbits <- nbits-1
     else nbits <- nbits+2
    spouse <- get.spouse(ped=ped,ID=id)
    children <- get.children(ped=ped,ID=id)
    siblings <- get.siblings(ped=ped,ID=id)
    parents <- get.parents(ped=ped,ID=id)
    print( ids <- unique(na.omit( c(spouse$ID, children$ID, siblings$ID, parents$ID)) ) )
    ids 
}

child.pedigrees <- function(ped, k=2) {
    children.ids <- ped[as.numeric(ped$children)==0,'ID']
    children.ids <- sample(children.ids,k)
    ids <- c()
    for (child.id in children.ids) {
    ids <- c(ids,child.id)
    parents <- get.parents(ID=child.id,ped=ped)
    ids <- c(ids,parents[,'ID'])
    while (nrow(parents)>0) {
    p <- parents[which(!as.logical(parents[,'founder'])),'ID']
    if (length(p)==0) break
    parents <- get.parents(ID=p,ped=ped)
    ids <- c(ids,parents[,'ID'])
    }
    }
    return(unique(ids))
}

ped.bits <- function(ped) {
     return( -length(which(as.logical(ped$founder))) + 2*length(which(!as.logical(ped$founder))) )
}

# Simple Mendelian rules.
# If parents are hom then child is hom.
# If parents are wt then child is wt.
complete.genotype.matrix <- function(ped,genotypes) {
    colnames(genotypes) <- gsub('.*_','',colnames(genotypes))
    print(samples <- intersect(ped$ID,colnames(genotypes)))
    geno <- genotypes[,samples]
    old <- rowSums(!is.na(geno))
    for (i in 1:nrow(geno)) {
        for (j in 1:ncol(geno)) {
            if (is.na(geno[i,j])) {
                father.id <- get.father(samples[j],ped)[,'ID']
                mother.id <- get.mother(samples[j],ped)[,'ID']
                if (is.na(father.id) | is.na(mother.id)) next
                dad.geno <- geno[i,father.id]
                mom.geno <- geno[i,mother.id]
                if (length(dad.geno)==0 | length(mom.geno)==0) next
                if (is.na(dad.geno) | is.na(mom.geno)) next
                if (dad.geno==0 & mom.geno==0) {
                    cat(i,samples[j],father.id,mother.id,':',0,'\n')
                    geno[i,j] <- 0 
                } else if (dad.geno==2 & mom.geno==2) {
                    cat(i,samples[j],':',2,'\n')
                    geno[i,j] <- 2 
               }
            }
        }
    }
  new <- rowSums(!is.na(geno))
  print(new-old)
  return(geno)
}


# all cases have the variant, closest founder parent has the variant
# founders are assumed to bring in variant to maximise penetrance
dominant.variant <- function(ped) {
    ped[,'dominant.variant'] <- 0
    # founders of affected must bring in the variant
    for (id in ped[which(ped$Affection==2),'ID']) {
        ped[id,'dominant.variant']<-1
        # if is.founder2 then continue
        if (is.founder2(id,ped)) next
        # if one of the parents already has the variant then continue
        p <- get.parents(id,ped)
        if (sum(ped[p$ID,'dominant.variant'])>0) next
        ped[get.founder(id,ped)[,'ID'],'dominant.variant']<-1
    }
    # remove ones for which we do not have DNA
    ped[which(is.na(ped$DNA)),'dominant.variant'] <- NA
    penetrance <- length(which(ped$Affection==2&!is.na(ped$DNA)))/length(which(ped$dominant.variant==1))
    cat('penetrance:',penetrance, '\n')
    return(ped)
}


# all cases have the variant, then variant is passed up to non-founders as it is assumed to be rare
# outside of this family.
# non-founders of affected must bring in the variant
dominant.rare.variant <- function(ped) {
    ped[,'dominant.variant'] <- 0
    for (id in ped[which(ped$Affection==2),'ID']) {
        ped[id,'dominant.variant']<-1
        # variant needs to propagate up to root
        while(TRUE) {
        # if is.founder2 then break
        if (is.founder2(id,ped)) break
        p <- get.parents(id,ped)
        # if both parents are founders then break
        # this should only happen once we reach the root of the pedigree
        if (is.founder2(p[1,'ID'],ped) & is.founder2(p[2,'ID'],ped)) break
        # assign variant to non-founder parent
        parent.id <- get.non.founder(id,ped)[,'ID']
        ped[parent.id,'dominant.variant']<-1
        # set parent.id as new id
        id <- parent.id
        }
    }
    # remove ones for which we do not have DNA
    ped[which(is.na(ped$DNA)),'dominant.variant'] <- NA
    penetrance <- length(which(ped$Affection==2&!is.na(ped$DNA)))/length(which(ped$dominant.variant==1))
    cat('penetrance:',penetrance, '\n')
    return(ped)
}




