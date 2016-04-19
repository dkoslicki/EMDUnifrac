import math,string,fpformat,random,re,sys,subprocess,getopt  # import of standard modules
import numpy as np
import dendropy

#print("WARNING: all branch lengths must be non-zero")
#print("WARNING: development still in progress. Use at your own risk")
#When doing all pairwise comparisons, I can get quite the speedup by making a single distance matrix from the tree, and just selecting
#the appropriate rows and columns based on the particular samples (the distance matrix creation is by far the most time consuming step).

class RandomTree:
	def __init__(self,alltips=10,nobr=0,pm=0.03,shape=0.5,mean=1):
		self.alltips=alltips   # number of tips
		self.nobr=nobr     # use branch lengths
		self.pm=pm         # probability of change per unit time
		self.shape=shape   # gamma shape parameter
		self.mean=mean     # mean of gamma dinstribution
	
	def constant_tree(self):  # function to generate a clock-like tree
		if self.alltips <=2:
			raise ValueError("error: need at least 3 tips.")
		tips=[]
		for i in range(1, self.alltips+1):
			tips.append("T"+str(i))
		Lb=[]
		for i in range(len(tips)):
			Lb.append(0)
		n=1
		dictionary={}
		while len(tips)!=1:
			R=random.random()
			tyme=(-(math.log(R))/len(tips))*self.pm
			fixtyme=fpformat.fix(tyme,5)
			brlens=float(fixtyme)
			for i in range(len(tips)):
				Lb[i]=Lb[i]+brlens
			nodeName = '@node%04i@' % n
			s1=random.choice(tips)
			i1=str(Lb[tips.index(s1)])
			del Lb[tips.index(s1)]
			tips.remove(s1)
			s2=random.choice(tips)
			i2=str(Lb[tips.index(s2)])
			del Lb[tips.index(s2)]
			tips.remove(s2)
			if self.nobr:
				nodo="("+s1+","+s2+")"
			else:
				nodo="("+s1+":"+i1+","+s2+":"+i2+")"
			dictionary[nodeName]=nodo
			tips.append(nodeName)
			Lb.append(0)
			n+=1
		findNodes=re.compile(r"@node.*?@", re.I) #to identify a node name
		lastNode = max(dictionary.keys())
		treestring = lastNode
		while 1:
			nodeList = findNodes.findall(treestring)
			if nodeList == []: break
			for element in nodeList:
				treestring=treestring.replace(element, dictionary[element])
		return treestring + ';'
	
	def variable_tree(self):  # function to generate a variable tree
		treestring=self.constant_tree()
		findbr=re.compile(":[0-9]+.[0-9]+[\),]")
		allbr=findbr.findall(treestring)
		dicbr={}
		for i in allbr:
			br=(i.split(':'))[1]
			brval=eval(br.strip('),'))
			beta=float(self.shape)/self.mean
			gammafactor=random.gammavariate(self.shape,beta)
			newbr=brval*gammafactor
			newbr1=fpformat.fix(newbr,5)
			dicbr[i]=newbr1
		for j in dicbr:
			if ',' in j:
				treestring=treestring.replace(j,':'+dicbr[j]+',') #This almost named all the internal nodes, but not quite. Try get non terminals
			elif ')' in j:
				treestring=treestring.replace(i,':'+dicbr[i]+')')
		return treestring


#Simulate environments, suitable for feeding directly to FastUnifrac. Will return distribution only on labeled nodes. Returns (envs)
def simulate_data(basis):
	#dtree = dendropy.Tree.get_from_string(tree_str, schema="newick", suppress_internal_node_taxa=False,store_tree_weights=True)
	#basis = [item.taxon.label for item in dtree.leaf_nodes()]
	#Use randomly generated samples on the leaves
	#Can change the distribution later if I want
	#cutoff = .2
	#Make sure we get a meaningful distribution
	#good_flag = 0
	#while good_flag == 0:
	weights_sample1 = np.random.exponential(scale=1.0, size=(1,len(basis))).transpose()
	#make it sparse
	#for i in range(len(weights_sample1)):
	#	if weights_sample1[i]<cutoff:
	#		weights_sample1[i] = 0
	#Multiply since fast_unifrac expects numbers > 1
	weights_sample1 = 1000*weights_sample1 
	weights_sample2 = np.random.exponential(scale=1.0, size=(1,len(basis))).transpose()
	#make it sparse
	#for i in range(len(weights_sample2)):
	#	if weights_sample2[i]<cutoff:
	#		weights_sample2[i] = 0
	weights_sample2 = 1000*weights_sample2
	#if weights_sample1.sum()>0 and weights_sample2.sum()>0:
	#	good_flag = 1
	envs = dict()
	i=0
	for node in basis:
		envs[node] = {'sample1':weights_sample1[i],'sample2':weights_sample2[i]}
		i = i+1
	return (envs)


def parse_tree(tree_str):
	# This function will parse a newick tree string and return the dictionary of ancestors T
	# and will also return the branch length dictionary l
	# NOTE: can probably clean this up so I compute Tint and lint directly instead of going over the tree twice
	#dtree = dendropy.Tree.get_from_string(tree_str, schema="newick", suppress_internal_node_taxa=False,store_tree_weights=True)
	dtree = dendropy.Tree.get(data=tree_str, schema="newick", suppress_internal_node_taxa=False,store_tree_weights=True)
	#Name all the internal nodes
	nodes = dtree.nodes()
	i=0
	for node in nodes:
		if node.taxon == None:
			node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
			i = i+1
	full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]  # i in path from root to j only if i>j
	full_nodes_in_order.reverse()
	nodes_in_order = [item.taxon.label for item in full_nodes_in_order]  # i in path from root to j only if i>j
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = full_nodes_in_order[i]
		parent = node.parent_node
		if parent != None:
			Tint[i] = nodes_to_index[parent.taxon.label]
			lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
	return (Tint,lint,nodes_in_order)



def parse_tree_file(tree_str_file):
	# This function will parse a newick tree string and return the dictionary of ancestors T
	# and will also return the branch length dictionary l
	dtree = dendropy.Tree.get(path = tree_str_file, schema="newick", suppress_internal_node_taxa=False,store_tree_weights=True)
	#Name all the internal nodes
	nodes = dtree.nodes()
	i=0
	for node in nodes:
		if node.taxon == None:
			node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
			i = i+1
	full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]  # i in path from root to j only if i>j
	full_nodes_in_order.reverse()
	nodes_in_order = [item.taxon.label for item in full_nodes_in_order]  # i in path from root to j only if i>j
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = full_nodes_in_order[i]
		parent = node.parent_node
		if parent != None:
			Tint[i] = nodes_to_index[parent.taxon.label]
			lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
	return (Tint,lint,nodes_in_order)

def parse_envs(envs, nodes_in_order):
	#This function will parse the environments and return a dictionary of vectors on the tree basis (nodes_in_order)
	nodes_in_order_dict = dict(zip(nodes_in_order,range(len(nodes_in_order))))
	for node in envs.keys():
		if node not in nodes_in_order_dict:
			print("Warning: environments contain taxa " + node + " not present in given taxonomic tree. Ignoring")
	envs_prob_dict = dict()
	for i in range(len(nodes_in_order)):
		node = nodes_in_order[i]
		if node in envs:
			samples = envs[node].keys()
			for sample in samples:
				if sample not in envs_prob_dict:
					envs_prob_dict[sample] = np.zeros(len(nodes_in_order))
					envs_prob_dict[sample][i] = envs[node][sample]
				else:
					envs_prob_dict[sample][i] = envs[node][sample]
	#Normalize samples
	samples = envs_prob_dict.keys()
	for sample in samples:
		envs_prob_dict[sample] = envs_prob_dict[sample]/envs_prob_dict[sample].sum()
	return (envs_prob_dict, samples)



#This will return the EMDUnifrac distance along with the flow
def EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q):
	num_nodes = len(nodes_in_order)
	F = dict()
	G = dict()
	Z = 0
	w = np.zeros(num_nodes)
	pos = dict()
	neg = dict()
	for i in range(num_nodes):
		pos[i] = set([])
		neg[i] = set([])
	for i in range(num_nodes):
		if P[i]>0 and Q[i]>0:
			F[(i,i)] = np.minimum(P[i],Q[i])
		G[(i,i)] = P[i] - Q[i]
		if P[i] > Q[i]:
			pos[i].add(i)
		elif P[i] < Q[i]:
			neg[i].add(i)
		posremove = set()
		negremove = set()
		for j in pos[i]:
			for k in neg[i]:
				if (j not in posremove) and (k not in negremove):
					val = np.minimum(G[(i,j)], -G[(i,k)])
					if val > 0:
						F[(j,k)] = np.minimum(G[(i,j)], -G[(i,k)])
						G[(i,j)] = G[(i,j)] - val
						G[(i,k)] = G[(i,k)] + val
						Z = Z + (w[j] + w[k])*val
					if G[(i,j)] == 0:
						posremove.add(j)
					if G[(i,k)] == 0:
						negremove.add(k)
		pos[i].difference_update(posremove)
		neg[i].difference_update(negremove)
		if i < num_nodes-1:
			for j in pos[i].union(neg[i]):
					if (Tint[i],j) in G:
						G[(Tint[i],j)] = G[(Tint[i],j)] + G[(i,j)]
					else:
						G[(Tint[i],j)] = G[(i,j)]
					w[j] = w[j] + lint[i,Tint[i]]
			pos[Tint[i]] |= pos[i]
			neg[Tint[i]] |= neg[i]
	return (Z, F)  # The returned flow is on the basis nodes_in_order and is given in sparse matrix dictionary format. eg {(0,0):.5,(1,2):.5}


#This will return the EMDUnifrac distance only
def EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q):
	num_nodes = len(nodes_in_order)
	Z = 0
	partial_sums = P - Q
	for i in range(num_nodes - 1):
		val = partial_sums[i]
		partial_sums[Tint[i]] += val
		Z += lint[i, Tint[i]]*abs(val)
	return Z




