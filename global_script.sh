#Structure of the code
for(i in 1:1000) {
	for(i in new_species_names) {
		if(species is already present in backbone tree) {do nothing}
		else {
			if(genus of new species is present) 
				{add new species at genus level}
			else {
				if (subfamily of new species is present) 
				{add new species at subfamily level}
				else {add at family level}
 }	} 		}	}

#WHAT SHOULD BE PREPARED BEFORE THE LOOP:
	#1. Get all needed data 
		#Data frame with species that we are aimimg to add and their taxonomy (prepared beforehand)
		species_to_add_taxonomy <- read.csv("./Data/species_to_add_taxonomy.csv", row.names=NULL)
		#Web tree taxonomy
		backbone_taxonomy <- read.csv("./Data/web_backbone_taxonomy.csv", row.names=NULL)
		#Upload the backbone tree
		backbone_tree <- read.tree(file = "./actinopt_12k_treePL.tre.xz")
	#2. Call libraries()
		library(ape)
		library(phytools)
		library(dplyr)
	#3. Create functions that we'll use in the loop
		#3.1 Function to obtain species taht belong to the certain genus of interest in the backbone tree
		get_species_from_genus <- function(tree, genus) {
			tips <- tree$tip.label # Get the tips (species) from the tree
			nodes <- which(sapply(strsplit(tips, "_"), `[`, 1) == genus) # Find the nodes that correspond to the given genus
			# Get the tips that are descendants of the nodes corresponding to the given genus
			  species <- character(0)
			  	for (node in nodes) {
				    descendants <- tree$tip.label[getDescendants(tree, node)]
				    species <- union(species, descendants)
			  		}
			return(species)
			}
		#3.2 Function to add species at the family level (halfway location)
		add_to_family <- function(tree, family, new_tip_label) {
	        mrca_node <- findMRCA(tree, family) #get mrca node
     	        daughter_nodes <- tree$edge[tree$edge[,1] == mrca_node, 2] # find descendant ('daughter) nodes from mrca
     		        insert_node <- sample(daughter_nodes, 1) #pick 1 daughter node
     		        edge_index <- which(tree$edge[,1] == mrca_node & tree$edge[,2] == insert_node) #find what is the index numb in the edge matrix that connects mrca & daughter node
     	        insert_position <- tree$edge.length[edge_index] / 2 #get insert position (halfway through the branch between mrca and daughter node)
     		        root_to_mrca <- node.depth.edgelength(tree)[mrca_node] #get distance from root ro mrca
     		        tree_length <- max(node.depth.edgelength(tree)) # get the length of the tree (from the root to the tip)
     	        edge_length_new <- tree_length - root_to_mrca - insert_position #get the edge length for the new species
     	        #add new species with 'bind.tip()'
     	        tree <- bind.tip(tree, new_tip_label, 
     	                              edge.length = edge_length_new, 
     	                              where = insert_node, 
     	                              position = insert_position) 
             return(tree)
           	} 
           
#GRANDE LOOP to get 1000 results of PD
PD_values <- numeric(100) #Create a numeric directory where the values of PD will be stored
new_species_names <- species_to_add_taxonomy$Species #get the list of Giovanni species

for(iteration in 1:100) {
	backbone_tree_1 <- backbone_tree #It is important to renew the backbone tree after each iteration (to get rid of previous additions)
	backbone_taxonomy <- backbone_taxon_SAVE #Also renew data frame with taxonomy
	#Loopy loop to add species to the tree
	for(i in new_species_names) {
		genus_n <- sub(" .*", "", i) #get the first word (genus) from the species name
		#Check if the species is present in backbone tree, if not, see 'else{}'
		if(i %in% backbone_taxonomy$Species) {
			#print(paste("Is already present in backbone tree:", i))
		}
		else {
			#Check if genus of new species is in the tree, if yes - add at genus level, if not, see 'else{}'
			if(genus_n %in% backbone_taxonomy$Genus) {
				backbone_genus <- get_species_from_genus(backbone_tree_1, genus_n)
					random_genus <- sample(backbone_genus, 1)
					edge_length_new <- backbone_tree_1$edge.length[which(backbone_tree_1$edge[,2] == which(backbone_tree_1$tip.label == random_genus))]
					add_name_1 <- gsub(" ", "_", i)
					backbone_tree_1 <- bind.tip(backbone_tree_1, add_name_1, 
				         						edge.length = edge_length_new/2, 
				         						where = which(backbone_tree_1$tip.label == random_genus), 
				         						position = edge_length_new / 2) 
				    #Important to bind taxonomic information about the new species to the backbone taxonomy to ensure a good downstream processing
				    backbone_taxonomy <- rbind(backbone_taxonomy, species_to_add_taxonomy[which(species_to_add_taxonomy$Species == i), ])
					#print(paste("New species successfully inserted within genus:", i)) 
	 				}
			else { 
				#get the subfamily information for new species
				new_subfamily <- species_to_add_taxonomy$Subfamily[which(species_to_add_taxonomy$Species == i)]
				#Check if this subfamily has representatives in backbone tree, if yes - bind at subfamily level, if not - see 'else{}' 
				if(new_subfamily %in% backbone_taxonomy$Subfamily) {
					subfamily.n <- backbone_taxonomy$Species[which(backbone_taxonomy$Subfamily == new_subfamily)] #get species from backbone tree that belong to the subfamily of interest
					subfamily_n <- gsub(" ", "_", subfamily.n) #ensure their name id written with '_' in between
					add_name_2 <- gsub(" ", "_", i) #ensure the name of new species is written with '_' in between
					backbone_tree_1 <- add_to_family(backbone_tree_1, subfamily_n, add_name_2)
					backbone_taxonomy <- rbind(backbone_taxonomy, species_to_add_taxonomy[which(species_to_add_taxonomy$Species == i), ])
					#print(paste("New species successfully inserted within subfamily:", i)) 
					}	
				else {
					#If the species was not bound neither at genus nor at subfamily level - it binds at family level
					new_family <- species_to_add_taxonomy$Family[which(species_to_add_taxonomy$Species == i)]
					family.n <- backbone_taxonomy$Species[which(backbone_taxonomy$Family == new_family)]
					family_n <- gsub(" ", "_", family.n)
					add_name_3 <- gsub(" ", "_", i)
					backbone_tree_1 <- add_to_family(backbone_tree_1, family_n, add_name_3)
					backbone_taxonomy <- rbind(backbone_taxonomy, species_to_add_taxonomy[which(species_to_add_taxonomy$Species == i), ])
					#print(paste("New species successfully inserted within family:", i)) 
				} } } }
	#Calculate Phylogenetic Diversity for each iteration
	PD_values[iteration] <- sum(backbone_tree_1$edge.length)
	print(paste("PD is calculated in iteration:", iteration))
	}

#AFTER PD IS CALCULATED
	#1. Get the 'summary()' of PD_values 
	summary(PD_values)
	#2. Plot PD values versus iteration number
	ggplot(data = data.frame(iteration = 1:100, PD = PD_values), 
       aes(x = iteration, y = PD)) +
     geom_point() + 
     geom_smooth(method = "loess", se = TRUE) +
     labs(x = "Iteration", y = "Phylogenetic Diversity",
          title = "Phylogenetic Diversity Across Iterations") +
     theme_minimal()
