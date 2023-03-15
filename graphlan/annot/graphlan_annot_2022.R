## Creates the annotation file for graphlan 

###########################
###   IO
###########################

eggnogg.d = "../../fivepseq_projects/microbiome/data/graphlan/eggnog_2022/"
annot.d = "../../fivepseq_projects/microbiome/data/graphlan/annot"
annot.f = "../../fivepseq_projects/microbiome/data/graphlan/annot/myannotR_2022.txt"
log.f = "../../fivepseq_projects/microbiome/data/graphlan/annot/myannotR_2022.log"

annotate_genus = F

my.species.f = file.path(annot.d, "my.species_2022.txt")
my.genomes.all.f = file.path(annot.d, "my.genomes.all.txt")
taxon.stats.f = file.path(annot.d, "taxon.stats.txt")
mup.stats.f = file.path(annot.d, "mup.stats.txt")

MIN_LIBSIZE <<- 10^20
MAX_LIBSIZE <<- 0 


my.species = read.table(my.species.f, header = T, sep = "\t")
my.genomes.all = read.table(my.genomes.all.f, header = T, sep = "\t", fill = T)
taxon.stats = read.table(taxon.stats.f, header = T, sep = "\t", comment.char = "", check.names = F)
mup.stats = read.table(mup.stats.f, header = T, sep = "\t", comment.char = "", check.names = F, row.names = 1)

out.con = file(annot.f, "wb")
log.con = file(log.f, "w")

###########################
###   Functions
###########################

out = function(vec, first.line = F){
  if(first.line)
    write(paste(as.character(vec), collapse = "\t"), out.con, append = F)
  else
    write(paste(as.character(vec), collapse = "\t"), out.con, append = T)
}

my.log = function(vec, first.line = F){
  if(first.line)
    write(paste(as.character(vec), collapse = " "), log.con, append = F)
  else
    write(paste(as.character(vec), collapse = " "), log.con, append = T)
}

get.max.ind = function(g, taxon.stats){
  ls = taxon.stats$libsize[which(taxon.stats$Genus == g)]
  max.ind = which(taxon.stats$Genus == g)[which.max(ls)]
  
  return(max.ind)
}

get.max.mup.ind = function(g, taxon.stats){
  k1 = which(colnames(taxon.stats) == "mup-13")
  k2 = which(colnames(taxon.stats) == "mup-17")
  mup.mat = taxon.stats[which(taxon.stats$Genus == g), k1:k2]
  max.ind = which(taxon.stats$Genus == g)[which.max(apply(mup.mat, 1, max))]
  
  return(max.ind)
}

clade.marker.size = function(g, taxon.stats){
  if(!(g %in% taxon.stats$Genus))
    cat("No genus ", g, "\n")
  ls = taxon.stats$libsize[which(taxon.stats$Genus == g)]
  mls = max(ls)
  if(mls > MAX_LIBSIZE)
    MAX_LIBSIZE <<- max(ls)
  if(mls < MIN_LIBSIZE)
    MIN_LIBSIZE <<- max(ls)
  s = max(ls)/5000
  if(s > 200)
    s = 200
  
  return(s)
}

clean.name = function(name){
  name = sub(".", "", name, fixed = T)
  name = sub("[", "", name, fixed = T)
  name = sub("]", "", name, fixed = T)
  name = sub("unclassified ", "", name, fixed = T)
  name = sub("Caulobacter crescentus", "Caulobacter vibrioides", name)
  return(name)
}

species.found = function(s, targets){
  s.clean = clean.name(s)
  
  for(t in targets){
    t.clean = clean.name(t)
    if(grepl(s.clean, t.clean))
      return(T)
  }
  
  return(F)
}

###########################
###   General options
###########################

out(c("clade_marker_edge_width", 0.1), first.line = T)
out(c("clade_marker_size", 0))
out(c("branch_color", "#D2D2D2"))


###########################
###   Rings
###########################
rings.f = file.path(annot.d, "rings_2022.txt")
rings = read.table(rings.f, sep = "\t", header = T, comment.char = "")
apply(rings, 1, function(x){
  out(c("ring_color", x[1], x[3]))
  out(c("ring_height", x[1], x[4]))
})


###########################
###   Phylum annotations
###########################
phyla.f = file.path(annot.d, "phyla.txt")
phyla = read.table(phyla.f, sep = "\t", header = T, comment.char = "")
apply(phyla, 1, out)

###########################
###   Phylum colors
###########################

phyla.colors.f = file.path(annot.d, "phylum.colors.txt")
phyla.colors = read.table(phyla.colors.f, sep = "\t", header = T, comment.char = "")
apply(phyla.colors, 1, out)

###########################
###   Genus annots
###########################
for(i in 1:nrow(my.species)){
  g = paste0("g_", my.species$Genus[i])
  cl = paste0("c_", my.species$Class[i])
  o = paste0("o_", my.species$Order[i])
  f = paste0("f_", my.species$Family[i])
  p = paste0("p_",my.species$Phylum[i])
  c = as.character(phyla.colors$color[which(phyla.colors$phylum == p)])
  
  genus = as.character(my.species$Genus[i])
  s = clade.marker.size(genus, taxon.stats)
  
  
  out(c(g, "annotation_background_color", c))
  if(annotate_genus)
    if(s > 20)
      out(c(g, "annotation", paste0(substr(genus, 1,7), ":", genus))) # uncheck this for clean image without annot
  
  out(c(g, "clade_marker_color", c))
  out(c(f, "clade_marker_color", c))
  out(c(o, "clade_marker_color", c))
  out(c(cl, "clade_marker_color", c))
  out(c(p, "clade_marker_color", c))
  out(c(g, "clade_marker_size", s))
  
}

###########################
###   Rings - periodicity
###########################


for(genus in unique(my.species$Genus)){
  
  g = paste0("g_", genus)
  max.ind = get.max.ind(as.character(genus), taxon.stats)
  
  # coverage
  j = which(rings$name == "coverage")
  ind = rings$index[j]
  h = rings$height[j]
  cov = log(taxon.stats$libsize[max.ind])/log(max(taxon.stats$libsize))
  out(c(g, "ring_alpha",  ind, h*cov))

  # periodicity 
  j = which(rings$name == "periodicity")
  ind = rings$index[j]
  h = rings$height[j]
  #out(c(g, "ring_height", ind, h*taxon.stats$periodicity[max.ind]))
  out(c(g, "ring_alpha", ind, h*taxon.stats$periodicity[max.ind]))
  
  # separator 1
  j = which(rings$name == "sep_white1")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
  j = which(rings$name == "sep_gray1")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
  j = which(rings$name == "sep_white2")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
  # frames
  j = which(rings$name == "F0")
  f0.ind = rings$index[j]
  h0 = rings$height[j]
  f0 = taxon.stats$F0[max.ind]
  
  
  j = which(rings$name == "F1")
  f1.ind = rings$index[j]
  h1 = rings$height[j]
  f1 = taxon.stats$F1[max.ind]
  
  j = which(rings$name == "F2")
  f2.ind = rings$index[j]
  h2 = rings$height[j]
  f2 = taxon.stats$F2[max.ind]
  
  if(taxon.stats$frame_pref[max.ind] == 0){
    lf0 = lf1 = lf2 = 0 
  } else if (taxon.stats$frame_pref[max.ind] < 0) {
    lf0 = lf1 = lf2 = 1
    k = which.min(c(f0,f1,f2))
    if(k == 1)
      lf0 = 0 
    else if(k == 2)
      lf1 = 0
    else 
      lf2 = 0
  } else {
    lf0 = lf1 = lf2 = 0
    k = which.max(c(f0,f1,f2))
    if(k == 1)
      lf0 = 1 
    else if(k == 2)
      lf1 = 1
    else 
      lf2 = 1
  }
  
  out(c(g, "ring_height", f0.ind, h0*lf0))
  out(c(g, "ring_height", f1.ind, h1*lf1))
  out(c(g, "ring_height", f2.ind, h2*lf2))
  
  # separator 2
  
  j = which(rings$name == "sep_white3")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
  j = which(rings$name == "sep_gray2")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
  j = which(rings$name == "sep_white4")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
}

###########################
###   Rings - enzymes
###########################
my.log(c("Processing enzymes"), first.line = T)

last.ind = ind + 1 
num.enzymes = 3
for(ind in last.ind : (last.ind + num.enzymes - 1)){
  enzyme = rings$name[ind]
  
  my.log(c("Enzyme: ", as.character(enzyme)))
  my.log(c("Color: ", as.character(rings$color[ind])))
  my.log(c("Ind: ", ind))
  
  eggnogg.f = file.path(eggnogg.d, enzyme, "extended_members.txt")
  eggnogg.df = read.table(eggnogg.f, sep = "\t", quote = "")
  
  # reduce genus names to the first part 
  inferred.genuses = unlist(lapply(eggnogg.df$V3, function(x){strsplit(as.character(x), split = " ")[[1]][1]}))
  my.log(c("Inferred eggnog genuses:", length(unique(inferred.genuses))))
  # count the fraction of appearance of that genus among my list of genuses 
  for(genus in unique(inferred.genuses)){
    i = which(inferred.genuses == genus)
    if(genus == "Synechocystis")
      genus = "Geminocystis"
    else if (genus == "Geminocystis")
      next
    
    if(genus %in% my.species$Genus){
      my.log(c("Genus: ", genus))
      g = paste0("g_", genus)
      n = 0
      t = sum(my.species$Genus == genus)
      j = which(my.species$Genus == genus)
      my.log("Species:")
      for(jj in j){
        s = as.character(my.species$Species[jj])
        if(species.found(s, as.character(eggnogg.df$V3[i]))){
          n = n + 1
          my.log(c("Found: ", s))
        } else {
          my.log(c("Not found: ", s))
        }
      }
      
      my.log(c("Total found:", n, "(out of ", t, ")"))
      
      alpha = n/t
      if(alpha > 1)
        alpha = 1
      out(c(g, "ring_height", ind, rings$height[ind]))
      out(c(g, "ring_alpha", ind, alpha))
    }
  }
  
  #  i = which(my.genomes.all$species %in% eggnogg.df$V3)
  #  genuses = unique(intersect(my.species$Genus, inferred.genuses))
  
  #  for(genus in genuses){
  #    g = paste0("g_", genus)
  #    out(c(g, "ring_height", ind, rings$height[ind]))
  #  }
}

###########################
###   Rings - mup
###########################
my.log(c("Processing mup pauses"), first.line = T)


for(genus in unique(my.species$Genus)){
  
  g = paste0("g_", genus)
  max.ind = get.max.mup.ind(as.character(genus), taxon.stats)
  
  # separator 
  j = which(rings$name == "sep_white5")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
  j = which(rings$name == "sep_gray3")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
  j = which(rings$name == "sep_white6")
  ind = rings$index[j]
  h = rings$height[j]
  out(c(g, "ring_height", ind, h))
  
  # mup
  mup.poss = as.character(seq(-13, -17, -1))  
  for(k in mup.poss){
    mup.pos = paste0("mup", k)
    j = which(rings$name == mup.pos)
    ind = rings$index[j]
    h = rings$height[j]  
    pause = taxon.stats[max.ind, mup.pos]
    out(c(g, "ring_alpha",  ind, h*pause))
  }
  
}

# SOMETHING MIGHT BE WRONG WITH THE ENZYMES - CHECK

close.connection(out.con)
close.connection(log.con)
