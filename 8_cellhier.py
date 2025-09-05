'''''''''''''''''
step1: cellular neighborhood
'''''''''''''''''
import sys
sys.path.append('/share/home/qlab/qlab_cdj/project/project_10x_ec/analysis/')
from cellhier.general import *
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import skimage
import pickle
import seaborn as sns
from cellhier.knn_graph_neighborhood import Neighborhoods
from sklearn.cluster import MiniBatchKMeans

os.chdir('/share/home/qlab/qlab_cdj/project/project_10x_ec/analysis')
cells = pd.read_csv('./codex_position_metadata.csv')
cells.head()
cells[['seurat_cluster']].value_counts()

ks = [10,20,30]
cluster_col = 'celltype'
exclude_cols = ['']
sum_cols = list(cells[cluster_col].unique())
cluster_cols = [a for a in sum_cols if a not in exclude_cols]
keep_cols = ['Exp', 'x', 'y', 'cell_id', 'celltype', 'seurat_cluster']

Neigh2 = Neighborhoods(cells=cells, ks=ks, cluster_col=cluster_col, sum_cols=sum_cols, 
                       keep_cols=keep_cols, reg='Exp', X='x', Y='y', add_dummies=True)
celltype_windows = Neigh2.k_windows()
BASE = cells[keep_cols]

for k in [10]:
  for n_clusters in [15]:
    print('normalizing', k)
    X_ = celltype_windows[k][cluster_cols].values
    #     X = transform_major_cells(X_,thresh = .2)
    X = X_/X_.sum(1, keepdims=True)
    X[np.isnan(X)] = 0
    print('clustering - ws:{} n_clusters: {}'.format(k, n_clusters))
    col_name = 'cluster_{}_ws{}_clusts{}'.format(cluster_col, k, n_clusters)
    assert col_name not in BASE.columns
    BASE[col_name] = MiniBatchKMeans(random_state=5, n_clusters=n_clusters).fit(X).labels_
BASE.head()

colors = [(0.00784313725490196, 0.24313725490196078, 1.0),
 (1.0, 0.48627450980392156, 0.0),
 (0.10196078431372549, 0.788235294117647, 0.2196078431372549),
 (0.9098039215686274, 0.0, 0.043137254901960784),
 (0.5450980392156862, 0.16862745098039217, 0.8862745098039215),
 (0.6235294117647059, 0.2823529411764706, 0.0),
 (0.9450980392156862, 0.2980392156862745, 0.7568627450980392),
 (0.6392156862745098, 0.6392156862745098, 0.6392156862745098),
 (1.0, 0.7686274509803922, 0.0),
 (0.0, 0.8431372549019608, 1.0),
 (0.0, 0.40790465, 0.16444444),
 (0.40392157, 0.0, 0.05098039),
 (0.0, 0.8431372549019608, 1.0),
 (0.0, 0.0, 0.0)]
colors = ["#abfff1","#FFFFB3","#ddd53e","#FB8072","#fba89f","#ec9ed2","#d561dd","#d0b33d",
          "lightgrey","#ff0004","#8e3af4","#77f538","#373bbf","#5893ba","#b1cdf1"]

catplot(BASE, hue='cluster_celltype_ws10_clusts15', X='x', Y='y', 
        palette={i: colors[i] for i in range(len(colors))}, figsize=15, size=4)

exp_names = cells['Exp'].unique()
cn_celltype = BASE[['cluster_celltype_ws10_clusts15', 'celltype']].value_counts()
cn_celltype[(0,)]
cn_celltype[(1,)]
cn_celltype[(2,)]
cn_celltype[(3,)]
cn_celltype[(4,)]
cn_celltype[(5,)]
cn_celltype[(6,)]
cn_celltype[(7,)]
cn_celltype[(8,)]
cn_celltype[(9,)]
cn_celltype[(10,)]
cn_celltype[(11,)]
cn_celltype[(12,)]
cn_celltype[(13,)]
cn_celltype[(14,)]

BASE['cluster_celltype_ws10_clusts15'] = BASE['cluster_celltype_ws10_clusts15'].astype('string')
BASE.to_csv('cns_celltype_ws10_c15.csv')

BASE['cns'] = "unknown"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '0'),'cns'] = "CN_Plasma"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '1'),'cns'] = "CN_TU_TUp_H"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '2'),'cns'] = "CN_TU"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '3'),'cns'] = "CN_Fibro"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '4'),'cns'] = "CN_Myeloid"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '5'),'cns'] = "CN_Endo"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '6'),'cns'] = "CN_T"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '7'),'cns'] = "CN_TU_TUp_H"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '8'),'cns'] = "CN_T"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '9'),'cns'] = "CN_B"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '10'),'cns'] = "CN_TUmsln"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '11'),'cns'] = "CN_TUp_H"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '12'),'cns'] = "CN_H"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '13'),'cns'] = "CN_Fibro_Endop"
BASE.loc[(BASE.cluster_celltype_ws10_clusts15 == '14'),'cns'] = "CN_Fibro"






'''''''''''''''''
step2: SC
'''''''''''''''''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

ks = [100]
cluster_col = 'cns'
exclude_cols = ['']
sum_cols = list(BASE[cluster_col].unique())
cluster_cols = [a for a in sum_cols if a not in exclude_cols]
keep_cols = ['Exp', 'x', 'y', 'cell_id', 'celltype', 'seurat_cluster', 'cluster_celltype_ws10_clusts15','cns']

Neigh2 = Neighborhoods(cells=BASE, ks=ks, cluster_col=cluster_col, sum_cols=sum_cols, 
                       keep_cols=keep_cols, reg='Exp', X='x', Y='y', add_dummies=True)
celltype_windows = Neigh2.k_windows()

#l = ['7','LZ', 'DZ', '5', '8', '1', '3','6', '4', '2', '0']
#l = ['7', '4', '1', '10', '14', '6', '13', '0', '12', '11', '9', '2', '8','5','3']
l = ['CN_TU',
 'CN_TU_TUp_H',
 'CN_B',
 'CN_Fibro_Endop',
 'CN_T',
 'CN_Endo',
 'CN_Plasma',
 'CN_Fibro',
 'CN_Myeloid',
 'CN_TUmsln',
 'CN_TUp_H',
 'CN_H']
lmap = {j:i for i,j in enumerate(l)}
palt = {str(i):j for i,j in enumerate(sns.color_palette('bright',15))}
palt = {str(i):j for i,j in enumerate(colors)}
palt = {'CN_TU':"#abfff1",
        'CN_TU_TUp_H':"#FFFFB3",
        'CN_B':"yellow",
        'CN_Fibro_Endop':"#FB8072",
        'CN_T':"#373bbf",
        'CN_Endo':"#ec9ed2",
        'CN_Plasma':"#d561dd",
        'CN_Fibro':"#d0b33d",
        'CN_Myeloid':"lightgrey",
        'CN_TUmsln':"#ff0004",
        'CN_TUp_H':"#8e3af4",
        'CN_H':"#77f538"}
sns.palplot(palt.values())

'''
this code samples the GC CNs and projects them into barycentric coordinates
'''
w = celltype_windows[100]
wgc = w.loc[w.loc[:,['CN_TU','CN_B','CN_Myeloid']].sum(axis=1)>9,:]
idx = wgc.index.values
x = wgc.loc[:,['CN_TU','CN_B','CN_Myeloid']]
proj = np.array([[0,0],[np.cos(np.pi/3),np.sin(np.pi/3)], [1,0]])
coords = np.dot(x/100,proj)

plt.figure(figsize=(7,7))
jit = .002
cols = [palt[a] for a in wgc['cns']]
plt.scatter(coords[:,0]+jit*np.random.randn(len(coords)),coords[:,1]+jit*np.random.randn(len(coords)),s = 2,alpha = .5, c = cols)
plt.axis('off')


'''
this is the code that finds the minimal combination of CNs
required to make up a threshold percentage of assignments in a window
combinations are stored as a sorted tuple
'''
def get_thresh_simps(x,thresh):
    sorts = np.argsort(-x, axis = 1)
    x_sorted = -np.sort(-x, axis = 1)
    cumsums = np.cumsum(x_sorted,axis = 1)
    thresh_simps = pd.Series([tuple(sorted(sorts[i,:(1+j)])) for i,j in enumerate(np.argmax(cumsums>thresh,axis = 1))])
    return thresh_simps

x = w.loc[:, sum_cols].values/100
simps = get_thresh_simps(x, .9)
simp_freqs = simps.value_counts(normalize=True)
simp_sums = np.cumsum(simp_freqs)

w['combination'] = [tuple(sum_cols[a] for a in s) for s in simps]
# save w
w.to_csv('windows100_celltype_ws10_c15_with_combinations.csv')


# this shows what proportion (y) of the total cells are assigned to the top x combinations
plt.figure(figsize=(20,5))
plt.plot(simp_sums.values)
plt.xticks(range(0,300,50),range(0,300,50),rotation = 90,fontsize = 10)


selected_cells = simps[simps.isin(simp_sums[simp_sums <= .99].index.values)]
catplot(w.loc[selected_cells.index.values], hue='cns', X='x', Y='y', palette=palt, figsize=10)

selected_cells = simps[simps.isin(simp_sums[simp_sums <= .9].index.values)]
catplot(w.loc[selected_cells.index.values], hue='cns', X='x', Y='y', palette=palt, figsize=10)

selected_cells = simps[simps.isin(simp_sums[simp_freqs >= .001].index.values)]
catplot(w.loc[selected_cells.index.values], hue='cns', X='x', Y='y', palette=palt, figsize=10)

g = nx.DiGraph()
thresh_cumulative = .99
thresh_freq = .001
# selected_simps = simp_sums[simp_sums<=thresh_cumulative].index.values
selected_simps = simp_freqs[simp_freqs >= thresh_freq].index.values

'''
this builds the graph for the CN combination map
'''
for e0 in selected_simps:
    for e1 in selected_simps:
        if (set(list(e0)) < set(list(e1))) and (len(e1) == len(e0)+1):
            g.add_edge(e0, e1)

tops = simp_freqs[simp_freqs >= thresh_freq].sort_values(ascending=False).index.values.tolist()[:20]


'''
this plots the CN combination map
'''
draw = g
pos = nx.drawing.nx_pydot.graphviz_layout(draw, prog='dot')

plt.figure(figsize=(40,20))
for n in draw.nodes():
    col = 'black'
    if len(draw.in_edges(n))<len(n):
        col = 'black'
    plt.scatter(pos[n][0],pos[n][1]-5, s = simp_freqs[n]*10000, c = col, zorder = -1)
    if n in tops:
        plt.text(pos[n][0],pos[n][1]-7, '*', fontsize = 25, color = 'white', ha = 'center', va = 'center',zorder = 20)
    delta = 8
    #plot_sim((pos[n][0]+delta, pos[n][1]+delta),n, scale = 20,s = 200,text = True,fontsize = 15)
    plt.scatter([pos[n][0]]*len(n),[pos[n][1]+delta*(i+1) for i in range(len(n))],c = [palt[l[i]] for i in n] ,marker = 's', zorder = 5,s = 400)
    
        
j = 0
for e0,e1 in draw.edges():
    weight = 0.2
    alpha = .3
    if len(draw.in_edges(e1))<len(e1):
        color = 'black'
        lw =1
        weight = 0.4
        
#     if (e0,e1) in set(draw.out_edges(tuple(sorted([lmap['3'],lmap['1']])))):
#         j+=1
#         print(j)
#         color = 'green'
#         weight = 2
#         alpha = 1
        
    if (lmap['CN_TU'] in e0) and (lmap['CN_Myeloid'] not in e0) and (lmap['CN_Myeloid'] in e1):
        color = 'green'
        weight = 2
        alpha = 1

    plt.plot([pos[e0][0], pos[e1][0]],[pos[e0][1], pos[e1][1]], color = color, linewidth = weight,alpha = alpha,zorder = -10)

plt.axis('off')
plt.savefig('CNM.pdf')
plt.show()





'''''''''
fig 3  assembly rules
'''''''''
import sys
from cellhier.general import *
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import skimage
import pickle
import seaborn as sns
import networkx as nx
from sklearn.neighbors import NearestNeighbors
from scipy.sparse.csgraph import connected_components
import itertools
from skimage.color import label2rgb
from scipy import ndimage as ndi

c = w 
#plot the tissues
catplot(c, hue = 'cns',X = 'x', Y = 'y', figsize = 10, palette = palt)

exps = ['EC4','EC5','EC7','EC10','EC11','EC12']
cns = sum_cols
all_cns = cns

'''
this code converts the cell CN assignments into an image,
then uses scipy label function to identify instances, which are connected components of each CN 
'''
imgs = {}
lab_imgs = {}
for exp in exps:#['spleen']:#exps:
    exp_cells = c.loc[c['Exp']==exp]

    downsamp = 50
    x = (exp_cells['x']/downsamp).astype(int)
    y = (exp_cells['y']/downsamp).astype(int)
    y = max(y)-y
    xdim = 1+max(x)
    ydim = 1+max(y)
    imgs[exp] = {}
    lab_imgs[exp] = {}

    for cn in cns:
        
        cnx = x[exp_cells['cns']==cn]
        cny = y[exp_cells['cns']==cn]
        img = np.zeros((ydim,xdim),dtype = int)
        np.add.at(img, (cny, cnx),1)
        imgs[exp][cn]=img

        plt.figure(figsize=(20,20))
        plt.subplot(1,2,1)
        #if cn in ['0','1','2','3','5','6','8']:
         #   thresh = 0
          #  dilate = False
        #else:
         #   thresh = 0
          #  dilate = True
        thresh = 0
        dilate = True
         
        plt.imshow(img>thresh,cmap = 'Greys_r')
        plt.axis('off')
        if dilate:
            dil = ndi.binary_dilation(img>thresh)
        else:
            dil = img>thresh
        lab_img, _ = ndi.label(dil)#ndi.binary_dilation(img))
        sizes = np.zeros((1+np.max(lab_img)))
        np.add.at(sizes,lab_img.flatten(),1)
        lab_img[np.isin(lab_img,np.arange(len(sizes))[sizes<20])] = 0
        
#         lab_img, _ = ndi.label(ndi.binary_dilation(lab_img>0))#ndi.binary_dilation(img))
#         sizes = np.zeros((1+np.max(lab_img)))
#         np.add.at(sizes,lab_img.flatten(),1)
#         lab_img[np.isin(lab_img,np.arange(len(sizes))[sizes<20])] = 0
        
        
        lab_imgs[exp][cn]=lab_img
        plt.title(cn)
        plt.subplot(1,2,2)
        plt.imshow(label2rgb(lab_img))
        plt.axis('off')
        plt.title(cn)
        plt.show()
        
    

'''
this code builds a graph for each tissue, vertices correspond to instances of a CN and edges correspond to adjacency
'''
gs = {}# nx.Graph()
for exp in exps:
    g = nx.Graph()
    for i,cn1 in enumerate(cns):
        for j,cn2 in enumerate(cns):
            if i<=j:
                continue
            lab1 = lab_imgs[exp][cn1]
            lab2 = lab_imgs[exp][cn2]
            intrx = np.zeros((np.max(lab1)+1, np.max(lab2)+1))
            np.add.at(intrx, (lab1.flatten(),lab2.flatten()),1)
            intrx[0,:] = -1
            intrx[:,0] = -1
            indexes = np.where(intrx>0)
            for e0,e1 in zip(*indexes):
                g.add_edge((exp,cn1,e0),(exp,cn2,e1))
    gs[exp] = g.copy()
        
    
    

'''
this code takes a motif and finds all the extensions, keeping track of which instances can be extended
the input info it needs includes graphs and dicts that it updates
g: tissue graph
tc: motif
tc_to_c : dict motifs -> instances
L_graph : graph of motif extensions
L_to_c : dict edges of L_graph -> instances of source of edge that extend 
tc_to_e : dict motifs -> experiments they're found it

a motif (typically referred to as tc in the code) is stored as a tuple (A,B,C) where:
A is a sorted tuple of sorted pairs of CNs representing the adjacencies in the motif
B is a sorted tuple of CNs representing the CNs in the motif
C is a sorted tuple of sorted pairs of CN representing the 'NOT' adjacencies
'''

def extend_motif(g,tc,tc_to_c, L_to_c, L_graph,tc_to_e):
    ty_edges, ty_nodes = tc
    samp = tc_to_c[tc]
    if len(samp) >5000:
        samp = list(samp)
        idd = np.arange(len(samp))
        np.random.shuffle(idd)
        samp = [samp[i] for i in idd[:5000]]
        
    for tm in samp:
        for tm_node in tm:
            for e0,e1 in g.edges(tm_node):
                e = e0
                if e ==tm_node:
                    e = e1
                exp, new_ty, new_inst = e
                
                new_edge = tuple(sorted([e0[1],e1[1]]))
                if  new_edge in ty_edges:
                    continue
                
                
                
                ty_node_set = set(ty_nodes)
                
                if new_ty in ty_node_set:
                    kv = {i:j for i,j in zip(ty_nodes,tm)}
                    if kv[new_ty]!=e:
                        continue
                        
                        
                    
                
                
                
                ty_node_set.add(new_ty)
                new_ty_nodes = tuple(sorted(ty_node_set))
                idx = {tt:nod for tt,nod in zip(ty_nodes,tm)}
                idx[new_ty] = e
                new_tm = tuple([idx[j] for j in new_ty_nodes])
                
                
                
                
                ty_edge_set = set(ty_edges)
                ty_edge_set.add(new_edge)
                
#                 for f0, f1 in g.edges(e):
#                     if (f0 in new_tm) and (f1 in new_tm):
#                         ty_edge_set.add(tuple(sorted((f0[1],f1[1]))))
                
                new_ty_edges = tuple(sorted(ty_edge_set))
                
                new_tc = (new_ty_edges,new_ty_nodes)
                tc_to_c.setdefault(new_tc,set()).add(new_tm)
                
                L_to_c.setdefault((tc,new_tc), set()).add(tm)
                tc_to_e.setdefault(new_tc,set()).add(tm[0][0])
                L_graph.add_edge(tc,new_tc)
                    
                
                
        
'''
this function takes a graph and uses extend motifs above to find all the motifs, their instances, 
as well as which motifs extend
'''
def build_motifs(g):
    tc_to_c = {}
    tc_to_e = {}

    for e0,e1 in g.edges():
        exp0,cn0,i0 = e0
        exp1,cn1,i1 = e1
        tc = ((tuple(sorted([cn0,cn1])),), tuple(sorted([cn0,cn1])))
        tm = tuple(sorted([e0,e1], key = lambda x: x[1]))
        tc_to_c.setdefault(tc, set()).add(tm)
        tc_to_e.setdefault(tc, set()).add(exp0)

    #take conserved motifs built from chains and adding
    L_to_c = {}
    L_graph = nx.DiGraph()
    for tc in tc_to_c:
        L_graph.add_node(tc)


    #expand them with add node motif    
    for i in range(15):
        leaves = [t for t in L_graph.nodes() if L_graph.out_degree(t)==0 and (len(t[1])<=4)]
        print(i,'---')
        print(len(leaves))
        for tc in leaves:
            extend_motif(g,tc,tc_to_c, L_to_c, L_graph,tc_to_e)
    return tc_to_c,L_to_c, L_graph




tc_to_cs = {}
L_to_cs = {}
L_graphs = {}
for exp in exps:
    tc_to_cs[exp],L_to_cs[exp],L_graphs[exp] = build_motifs(gs[exp])





def pb_inst(c1,c2,inst_set):
    instmap = {cn:i for i,cn in enumerate(c2[1])}
    return set([tuple([inst2[instmap[cn]] for cn in c1[1]]) for inst2 in inst_set])
    
    
    
    
    
'''
this is the code that builds the graph of minimal implications 
''' 

def added_piece(e0,e1):
    if len(e1)==3:
        l = set(list(e1[2])).intersection(set(list(e1[1])))
        if len(l)==2:
            ad = 'not=_'+str(e1[2])
        if len(l)==1:
            ad = str(list(l)[0]) + '_not_' + str(e1[2])
    if len(e1)==2:
        new_edge = list(set(list(e1[0]))- set(list(e0[0])))[0]
        if len(e1[1])==len(e0[1]):
            ad = '='+str(new_edge)
        else:
            new_node = list(set(list(e1[1]))-set(list(e0[1])))[0]
            ad_point = list(set(list(new_edge))-set([new_node]))[0]
            ad = str(ad_point)+'_'+str(new_edge)
    return ad


hom_graphs = {}
pruned_graphs = {}
mod_graphs = {}
added_pieces = {}
for exp in exps:
    L_graph = L_graphs[exp].copy()
    tc_to_c = tc_to_cs[exp].copy()
    L_to_c = L_to_cs[exp].copy()

    for node in [a for a in L_graph.nodes() if len(a[1])==2]:
        for i,x in enumerate(node[1]):
            tc_to_c[(tuple(),(x,))] = set()
            L_to_c[((tuple(),(x,)),node)] = set()


    for node in [a for a in L_graph.nodes() if len(a[1])==2]:
        for i,x in enumerate(node[1]):
            L_graph.add_edge((tuple(),(x,)),node)
            for ch in tc_to_c[node]:
                tc_to_c.setdefault((tuple(),(x,)), set()).add(ch[i])
                L_to_c.setdefault(((tuple(),(x,)),node), set()).add(ch[i])


    #this bit adds the 'nots' for all possible cn additions (although only for individual nodes)
    for e0 in list(L_graph.nodes()):
        if len(e0[0])>0:
            continue
        for n in e0[1]:
            for cn in all_cns:
                if cn==n:
                    continue
                if ((),(cn,)) not in tc_to_c:
                    continue
                
                added_comp = tuple(sorted([n,cn]))
                comp = (e0[0], e0[1], added_comp)
                L_graph.add_edge(e0,comp)
                new_ty = (tuple(sorted(set([added_comp] + list(e0[0])))), tuple(sorted(set([cn] + list(e0[1])))))
                if (e0,new_ty) in L_to_c:
                    tc_to_c[comp] = tc_to_c[e0] - L_to_c[e0,new_ty]
                    L_to_c[e0,comp] = tc_to_c[e0] - L_to_c[e0,new_ty]
                else:
                    tc_to_c[comp] = tc_to_c[e0]
                    L_to_c[e0,comp] = tc_to_c[e0]
            

            
            
    #this bit adds the 'nots' for only those additions at each step for which there are actual lifts
    for e0,e1 in list(L_graph.edges()):
        if len(e0[0])<1:
            continue
        added_edge = list(set(list(e1[0])) - set(list(e0[0])))[0]
        comp = (e0[0], e0[1], added_edge)
        L_graph.add_edge(e0,comp)
        tc_to_c[comp] = tc_to_c[e0]-L_to_c[e0,e1]
        L_to_c[e0,comp] = tc_to_c[e0]-L_to_c[e0,e1]
    
    mod_graphs[exp] = L_graph.copy()
    
    # pruned graph is the graph of all assembly rules
    pruned_graph = nx.DiGraph()
    thresh = .7
    for e in L_graph.edges():
        if len(tc_to_c[e[0]])<5:
            continue
            
        if len(e[1])==3:
            k2_0 = tuple(sorted(set(list(e[0][0])).union(set([e[1][2]]))))
            k2_1 = tuple(sorted(set(list(e[0][1])).union(set(list(e[1][2])))))
            k2 = (k2_0, k2_1)
            if len(tc_to_c.setdefault(k2, set())) >5:
                continue
        if (len(L_to_c[e])/len(tc_to_c[e[0]]))>thresh:
            pruned_graph.add_edge(e[0],e[1])
    
    #this bit deletes edges which aren't transitive
    del_edges = set()
    thresh2 = .7
    for e0,e1 in pruned_graph.edges():
        for e00 in nx.ancestors(pruned_graph,e0):
            if len(pb_inst(e00,e1,tc_to_c[e1]))/len(tc_to_c[e00]) < thresh2:
                del_edges.add((e0,e1))
    
    for e0,e1 in del_edges:
        pruned_graph.remove_edge(e0,e1)
    

    # hom graph is the graph of basic rules
    hom_graph = nx.DiGraph()            
    thresh3 = 0.3
    for i in range(10):
        nodes = [n for n in pruned_graph.nodes() if ((len(n)==3) and len(n[0])==(i-1)) or ((len(n)==2) and len(n[0])==i)]

        for node in nodes:
            for _,e1 in pruned_graph.out_edges(node):
                add = True
                diff = added_piece(node,e1)
                for d in nx.ancestors(L_graph,node):
                    if d not in pruned_graph.nodes():
                        continue
                    for _,w in pruned_graph.out_edges(d):
                        if added_piece(d,w)==diff:
                            add = False
                if add:
                    hom_graph.add_edge(node,e1)
                    
    print(len(hom_graph.edges()))
    pruned_graphs[exp] = pruned_graph.copy()
    hom_graphs[exp] = hom_graph.copy()



    added_pieces[exp] = {}
    for e in hom_graphs[exp].edges():
        added_pieces[exp].setdefault(added_piece(*e),set()).add(e)




'''
this just counts all the extensions and rules of each kind
'''
counts = {}
count_sets0 = {}
count_sets1 = {}
counts_all = {}
graphs = {}
for exp in exps:
    counts[exp] = {}
    count_sets0[exp] = {}
    count_sets1[exp] = {}
    counts_all[exp]= {}

    for e0, e1 in hom_graphs[exp].edges():
        if len(e1)==2:
            g0 = nx.Graph()
            for a in e0[1]:
                g0.add_node(a)
            for a0,a1 in e0[0]:
                g0.add_edge(a0,a1)
                
            g1 = nx.Graph()
            for a in e1[1]:
                g1.add_node(a)
            for a0,a1 in e1[0]:
                g1.add_edge(a0,a1)
                
            key = (tuple(sorted([g0.degree(n) for n in g0.nodes()])),tuple(sorted([g1.degree(n) for n in g1.nodes()])))
            graphs[key] = (g0,g1)
            counts[exp].setdefault(key,0)
            count_sets0[exp].setdefault(key,set()).add(e0)
            count_sets1[exp].setdefault(key,set()).add(e1)
            counts_all[exp].setdefault(key,set()).add((e0,e1))
            
            counts[exp][key]+=1




rule_types = set(counts['EC4'].keys()).union(set(counts['EC5'].keys())).union(set(counts['EC7'].keys())).union(set(counts['EC10'].keys())).union(set(counts['EC11'].keys())).union(set(counts['EC12'].keys()))
rule_types = sorted(rule_types,key = lambda x: (len(x[0]), x[0], len(x[1]),x[1]))





plt.figure(figsize=(21,1))
for j,key in enumerate(rule_types):
    plt.subplot(1,21,j+1)
    g0,g1 = graphs[key]
    pos = nx.drawing.nx_pydot.graphviz_layout(g1, prog='neato')
    for i,n in enumerate(g1.nodes()):
        col = sns.color_palette('tab10')[i]
        plt.scatter(pos[n][0],pos[n][1], s = 40, c = 'black', zorder = 10)
        #plt.text(pos[n][0], pos[n][1], ['A','B','C','D','E'][i])

    
    for e0,e1 in g0.edges():
        plt.plot([pos[e0][0], pos[e1][0]], [pos[e0][1], pos[e1][1]], color = 'black')
        
    e0,e1 = list(set(g1.edges)-set(g0.edges))[0]
    plt.plot([pos[e0][0], pos[e1][0]], [pos[e0][1], pos[e1][1]],linewidth = 5, color = 'green')

    
    plt.axis('off')
plt.show() 
    





l=0
plt.figure(figsize=(21,2))
for exp in exps:
    for _,key in enumerate(rule_types):
        l+=1
        plt.subplot(6,21,l)
        plt.text(0,0,  str(counts[exp].setdefault(key,0)),fontsize = 15, ha = 'center', va = 'center')# + ' / ' + str(len(count_sets0[exp].setdefault(key,set()))) + ' / ' + str(len(count_sets1[exp].setdefault(key,set()))))
        plt.axis('off')
    

    



'''
this code plots the graphs as shown in the figure, adding the edges where a rule is not present
'''
for exp in exps:
    hom_graph = hom_graphs[exp]
    sub = nx.DiGraph()
    for e in hom_graph.edges():
        if (len(e[0][1])<=1):# or ((len(e[0][1])==3) and (len(e[1][1])==3)):
            sub.add_edge(e[0],e[1], col = 'blue')
            if len(e[1][1])==2:
                n = list(set(list(e[1][1])) - set(list(e[0][1])))[0]
                sub.add_edge(e[1],((),(n,)), col = 'black')

    draw = sub    

    draw2 = nx.DiGraph()
    nd = {e1:i for i,e1 in enumerate(draw.nodes())}
    ndinv = {i:e1 for i,e1 in enumerate(draw.nodes())}
    for n in draw.nodes():
        draw2.add_node(nd[n])
    for e in draw.edges():
        draw2.add_edge(nd[e[0]],nd[e[1]])
    pos = nx.drawing.nx_pydot.graphviz_layout(draw2, prog='dot')
    pos = {ndinv[k]:v for k,v in pos.items()}
    print(len(pos.keys()))
    print(len(draw2.nodes()))
    print(len(draw.nodes()))





    plt.figure(figsize=(15,15))
    for n in draw.nodes():
        col = 'black'
        if len(n[1])==1:
            if len(n)==2:
                plt.scatter(pos[n][0], pos[n][1], s = 100, c = [palt[n[1][0]]])
            if len(n)==3:
                plt.scatter(pos[n][0]-10, pos[n][1], s = 100, c = [palt[n[1][0]]])
                ne = list(set(list(n[2])) - set([n[1][0]]))[0]
                plt.scatter(pos[n][0]+10, pos[n][1], s = 100, c = [palt[ne]])
                plt.plot([pos[n][0]-10,pos[n][0]+10], [pos[n][1],pos[n][1]], color= 'red')
                
        if len(n[1])==2:
            plt.scatter(pos[n][0]-10, pos[n][1], s = 100, c = [palt[n[0][0][0]]])
            plt.scatter(pos[n][0]+10, pos[n][1], s = 100, c = [palt[n[0][0][1]]])
            plt.plot([pos[n][0]-10,pos[n][0]+10], [pos[n][1],pos[n][1]], color= 'grey')
    cols = nx.get_edge_attributes(draw, 'col')        
    for e0,e1 in draw.edges():
        weight = 0.4
        zorder = 5
        col = cols[e0,e1]
        if col =='black':
            zorder = 4
        


        x = pos[e0][0]
        dx = pos[e1][0]-x
        if pos[e1][1]<pos[e0][1]:
            y = pos[e0][1]-20
            dy = pos[e1][1]+20-y
        else:
            y = pos[e0][1]+20
            dy = pos[e1][1]-20-y

        plt.arrow(x,y,dx,dy,head_length =5, head_width = 5,length_includes_head= True, color = col, zorder = zorder)
    plt.title(exp)
    plt.axis('off')
    plt.savefig(f"assembly_rules_{exp}.pdf")
    plt.show()










'''''''''
fig4
'''''''''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from tqdm.notebook import tqdm
from statsmodels.stats.multitest import fdrcorrection

#get windows for CRC dataset
cells2 = pd.read_pickle('submissiondata/cells2_salil')
m = cells2.groupby('neighborhood10').apply(lambda x: x['ClusterName'].value_counts(normalize = True)).unstack()

'''
this is the code that finds the minimal combination of CNs
required to make up a threshold percentage of assignments in a window
combinations are stored as a sorted tuple
'''
def get_thresh_simps(x,thresh):
    sorts = np.argsort(-x, axis = 1)
    x_sorted = -np.sort(-x, axis = 1)
    cumsums = np.cumsum(x_sorted,axis = 1)
    thresh_simps = pd.Series([tuple(sorted(sorts[i,:(1+j)])) for i,j in enumerate(np.argmax(cumsums>thresh,axis = 1))])
    return thresh_simps


w = pd.read_csv('submissiondata/windows_and_SCs_spotsover750cells.csv')
w['neighborhood10'] = w['neighborhood10'].astype('category')
l = ['0','1','2','3','4','5','6','7','8','9']
k = 50
x = w.loc[:,l].values/k
simps = get_thresh_simps(x,.90)
w['combination'] = [tuple(l[a] for a in s) for s in simps]
bulk_simp_freqs = (w.groupby('combination').size()/len(w)).to_dict()


spot_SCs = w.groupby(['patients','combination']).size().reset_index()
sc_pat_counts = spot_SCs.loc[spot_SCs[0]>25].loc[:,['patients','combination']].groupby('combination').size()
selected_simps = sc_pat_counts[sc_pat_counts>=10].index.values



checkpoint_cells = ['CD11c+Ki67+','CD11c+PD-1+','CD11c+VISTA+',
'CD68+CD163+PD-1+','CD68+CD163+Ki67+','CD68+CD163+VISTA+',
'CD38+Ki67+',
'CD20+PD-1+','CD20+Ki67+',
'CD8+PD-1+', 'CD8+Ki67+',
'Treg-PD-1+','Treg-Ki67+',
'CD4+PD-1+', 'CD4+Ki67+']

parents = ['CD4+', 'CD8+', 'CD20+', 'CD38+', 'CD68+CD163+', 'CD11c+', 'CD25+FOXP3+']


parent_checkpointlist = {'CD11c+':['CD11c+Ki67+','CD11c+PD-1+','CD11c+VISTA+'],
'CD68+CD163+': ['CD68+CD163+PD-1+','CD68+CD163+Ki67+','CD68+CD163+VISTA+'],
'CD38+': ['CD38+Ki67+'],
'CD20+':['CD20+PD-1+','CD20+Ki67+'],
'CD8+':['CD8+PD-1+', 'CD8+Ki67+'],
'CD25+FOXP3+':['Treg-PD-1+','Treg-Ki67+'],
'CD4+':['CD4+PD-1+', 'CD4+Ki67+']}

all_spots = w['spots'].unique()



parent_spot_cells = {}
for parent in parents:
    parent_spot_cells[parent] = w[w[parent]==1].groupby('spots').groups


parent_spot_perms = {}
num_perm = 30000

np.random.seed(4)
for parent in parents:
    parent_spot_perms[parent] = {}
    for spot in tqdm(all_spots):
        spot_parent_idxs = parent_spot_cells[parent][spot].values
        parent_spot_perms[parent][spot] = np.zeros((num_perm, len(spot_parent_idxs)),dtype = np.int64)
        for i in range(num_perm):
            shuf = spot_parent_idxs.copy()
            np.random.shuffle(shuf)
            parent_spot_perms[parent][spot][i] = shuf



ens = {}
pvals = {}
for sc in tqdm(selected_simps):
    ens[sc] = {}
    pvals[sc] = {}
    sc_cells = w[w['combination']==sc]
    for parent,checkpointlist in parent_checkpointlist.items():
        sc_parent_cells = sc_cells[sc_cells[parent]==1]
        sc_spot_parent_counts =sc_parent_cells.groupby('spots').size()
        ct_df = w.loc[:,checkpointlist]
        samp_array = np.zeros((len(sc_spot_parent_counts), num_perm, len(ct_df.columns)))
        for idx,spot in enumerate(sc_spot_parent_counts.index.values):
            if sc_spot_parent_counts.loc[spot] == parent_spot_perms[parent][spot].shape[1]:
                start ==0
            else:
                start = np.random.randint(parent_spot_perms[parent][spot].shape[1]-sc_spot_parent_counts.loc[spot])
                
            cells_to_get = parent_spot_perms[parent][spot][:,start:(start+sc_spot_parent_counts.loc[spot])]
            samp_array[idx] = ct_df.loc[cells_to_get.flatten()].values.reshape((cells_to_get.shape[0], cells_to_get.shape[1],len(ct_df.columns))).sum(axis = 1)
        samp_mean = samp_array.sum(axis = 0)/len(sc_parent_cells)
        obs_mean = ct_df.loc[sc_parent_cells.index.values].mean(axis = 0)
        p_greater = 1 - (samp_mean < obs_mean[None,:]).mean(axis = 0)
        p_less = 1 - (samp_mean > obs_mean[None,:]).mean(axis = 0)
        ens[sc,parent] = np.log2(((1e-4+obs_mean[None,:])/(1e-4+samp_mean)))
        pvals[sc,parent] = np.maximum(1/num_perm, np.minimum(p_greater,p_less))
        
        

to_correct = []
names = []
for sc in selected_simps:
    sc_cells = w[w['combination']==sc]
    for parent,checkpointlist in parent_checkpointlist.items():
        if parent =='CD38+':
            continue
        sc_parent_cells = sc_cells[sc_cells[parent]==1]
        sc_pat_parent_counts =sc_parent_cells.groupby('patients').size()
        if len(sc_pat_parent_counts)>=5:
            en = np.median(ens[sc,parent],axis = 0)
            for j,chk in enumerate(checkpointlist):
                if en[j]>0:
                    to_correct.append(pvals[sc,parent][j])
                    names.append((sc,chk,'pos'))
                if en[j]<0:
                    to_correct.append(pvals[sc,parent][j])
                    names.append((sc,chk,'neg'))




chk_par  = {'CD11c+Ki67+': 'CD11c+',
 'CD11c+PD-1+': 'CD11c+',
 'CD11c+VISTA+': 'CD11c+',
 'CD68+CD163+PD-1+': 'CD68+CD163+',
 'CD68+CD163+Ki67+': 'CD68+CD163+',
 'CD68+CD163+VISTA+': 'CD68+CD163+',
 'CD38+Ki67+': 'CD38+',
 'CD20+PD-1+': 'CD20+',
 'CD20+Ki67+': 'CD20+',
 'CD8+PD-1+': 'CD8+',
 'CD8+Ki67+': 'CD8+',
 'Treg-PD-1+': 'CD25+FOXP3+',
 'Treg-Ki67+': 'CD25+FOXP3+',
 'CD4+PD-1+': 'CD4+',
 'CD4+Ki67+': 'CD4+'}


parent_names = {
 'CD4+':'CD4$^+$ T: ', 
'CD8+': 'CD8$^+$ T: ',
    'CD20+': 'B cell: ',
    'CD38+': 'plasma/activated cell: ',
    'CD68+CD163+': 'MΦ: ',
    'CD11c+': 'DC: ',
    'CD25+FOXP3+': 'Treg: '   
}


sc_to_sig_chks = {}
for sc, chk, ud in [names[j] for j in np.where(fdrcorrection(to_correct)[1]<0.05)[0]]:
    sc_to_sig_chks.setdefault(sc,{}).setdefault(parent_names[chk_par[chk]], {}).setdefault(ud, set()).add(chk.split(chk_par[chk])[-1][:-1])
    


cn_name_map = {'0':1}



'''
this builds the graph for the CN combination map
'''
g = nx.DiGraph()
for e0 in selected_simps:
    for e1 in selected_simps:
        if (set(list(e0))<set(list(e1))) and (len(e1) == len(e0)+1):
            g.add_edge(e0,e1)
            
            
            
            
            

draw = g
pos = nx.drawing.nx_pydot.graphviz_layout(draw, prog='dot')
plt.figure(figsize = (60,30))
pal = sns.color_palette('bright',10)
for n in draw.nodes():
    colpal = sns.color_palette('bwr',100)
    vmin = -2
    vmax = 2
    
    col = 'black'
    if len(draw.in_edges(n))<len(n):
        col = 'black'
  
    plt.scatter(pos[n][0],pos[n][1]-5, s = 5*bulk_simp_freqs[n]*10000, c = ['black'], vmax = 20, zorder = -1)
    
    
    delta = 5

    plt.scatter([pos[n][0]]*len(n),[pos[n][1]+delta*(i+1) for i in range(len(n))],c = [pal[int(i)] for i in n] ,marker = 's', alpha = .3,zorder = 5,s = 400)
    for i in range(len(n)):
        plt.text(pos[n][0],pos[n][1]+delta*(i+1),cn_name_map.setdefault(n[i],n[i]), ha = 'center', va = 'center',fontsize = 18,zorder = 20,weight = 'bold')
    
    

    rrr = 1
    for ct in sorted(sc_to_sig_chks.setdefault(n, {}).keys()):
        plt.text(pos[n][0],-8-5*rrr+ pos[n][1],ct,fontsize = 20, ha = 'right', va = 'center', c = 'black')
        for mark in sorted(sc_to_sig_chks[n][ct].setdefault('pos', set())):
            plt.text(pos[n][0],-8-5*rrr+ pos[n][1],mark,fontsize = 20, ha = 'left', va = 'center', c = 'red')
            rrr +=1 
        for mark in sorted(sc_to_sig_chks[n][ct].setdefault('neg', set())):
            plt.text(pos[n][0],-8-5*rrr+ pos[n][1],mark,fontsize = 20, ha = 'left', va = 'center', c = 'blue')
            rrr +=1 
            

        
j = 0
for e0,e1 in draw.edges():
    weight = .8
    alpha = .3
    
    color = 'black'


    plt.plot([pos[e0][0], pos[e1][0]],[pos[e0][1], pos[e1][1]], color = color, linewidth = weight,alpha = alpha,zorder = -10)

plt.axis('off')

plt.show()



