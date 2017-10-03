import numpy as np
import matplotlib as mpl 
mpl.use('agg')
import matplotlib.pyplot as plt

def main(filename):
    #tossing open the first argument as the presumed csv
    import csv
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        header = reader.next()
        dataset=np.array([row for row in reader])
    
    from os.path import basename
    name=basename(filename)

    #is_ffpe should stay in 0th position
    col_names=['is_ffpe','insert_stdev','unmapped_reads','aligned_bases','reads_on_target','average_read_length','reads_per_start_point','insert_mean','total_reads','soft_clip_bases']
    col_eyes=map(lambda i: header.index(i), col_names)

    plotindex=1

    #make the figure letter-sized
    fig=plt.figure(figsize=(8.5,11))

    # 3 columns & figure out how many rows
    from math import ceil
    rows=ceil((len(col_names)-1)/3)
    cols=3

    #skip 'is_ffpe' and start at the next thing
    for i in range(1,len(col_eyes)):
        #current header value
        hindex=col_eyes[i]
        #slice out the is_ffpe and current header column
        tmp_arr=dataset[:,np.array([0,hindex])]

        #set up the subplot
        ax=fig.add_subplot(rows,cols,plotindex)
        plotindex+=1

        box_em_up(ax,header[hindex],tmp_arr.tolist())

    #fix the multiple plots so they don't overlap each other
    plt.tight_layout()
    fig.savefig('.'.join([name,'png']))


def box_em_up(ax,label,array):
    # pull out the rows that are ffpe and not ffpe (probably fresh frozen; may be LCM?)
    ffpe=filter(lambda x : x[0]=='True',array)
    ff=filter(lambda x : x[0]=='False',array)
    
    
    #get rid of the first column is_ffpe and leave only numbers
    ffpe_pts=np.array([numbers[1:] for numbers in ffpe]).astype(np.float)
    ff_pts=np.array([numbers[1:] for numbers in ff]).astype(np.float)

    ax.set_ylabel(label)
    bp=ax.boxplot([ffpe_pts,ff_pts],labels=['FFPE','Fresh Frozen'], patch_artist=True)
    make_boxes_colors(bp,'#ffffff')
    # test whether these two values are significant
    # if they are, plot them on the figure    
    from scipy.stats import ttest_ind
    pval=ttest_ind(ff_pts.tolist(),ffpe_pts.tolist()).pvalue
    if pval<0.01:
        labelval="=".join(['p','{:.2E}'.format(float(pval))])
        # list the p value
        ax.text(0.5,1,labelval, ha='center',va='top',transform=ax.transAxes)
        # add an asterisk
        ax.annotate(xycoords='data',xy=(1,np.amax(ffpe_pts)),s="*",weight='extra bold')
        make_boxes_colors(bp,'#1b9e77')

def make_boxes_colors(bp,color):
    for box in bp['boxes']:
        box.set( facecolor = color)


if __name__ == "__main__":
    import sys
    main(sys.argv[1])
