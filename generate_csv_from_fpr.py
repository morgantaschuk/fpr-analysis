#!/usr/bin/python

projects=['BART2', 'BartlettTS', 'CPC-GENE', 'CYT', 'DCIS', 'DYS', 'EACdysplasia', 'FFPER', 'GECCO', 'GECCOSequencing', 'GECCOTestSamples', 'GLCS', 'IDP', 'JDRT', 'JoshuaProstateSeq', 'Liposarcoma_Sequencing', 'MA5', 'MCNT', 'miRNASeq_Optimization', 'OCT', 'Ovarian_Brain_Colon_Exome_Seq', 'PCSI', 'SteveGallingerEarlyOnsetCRC', 'TGL01', 'TGL02', 'TGL03', 'TGL05', 'TGL07', 'TGL08', 'TGL09', 'THB', 'ValentinaLowInputTest']

def main(jsonfile):
    srli_libs, srli_workflows = read_fpr(jsonfile)
    bqcs=open_bamqcs(get_files(srli_workflows,'BamQC',set()))
    # first: insert size by ffpe and non-ffpe
    # second: % on target by ffpe and non-ffpe
    import collections
    gm_fields=['is_ffpe','multiplex',"target_size","aligned_bases", "insert_stdev", "soft_clip_bases", "reads_on_target", "average_read_length", "reads_per_start_point", "insert_mean", "total_reads", "unmapped_reads","run_name","lane","barcode"]
    GraphMe=collections.namedtuple('GraphMe', gm_fields)

    import csv
    import basename
    with open(".".join([basename(jsonfile),'csv']), 'wb') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=gm_fields)
        writer.writeheader()
        for b in bqcs:
            print(b)
            sr=b.run_name
            l=str(b.lane)
            i=b.barcode
            print(" ".join([sr,l,i]))
            if srli_libs[sr][l][i]==0:
               continue 
            i_f=srli_libs[sr][l][i].is_ffpe
            m=srli_libs[sr][l]['multiplex']
            gm=GraphMe(i_f,m,**b._asdict())
            writer.writerow(gm._asdict())



def open_bamqcs(bamqcs):
    fields=["target size","aligned bases", "insert stdev", "soft clip bases", "reads on target", "average read length", "reads per start point", "insert mean", "total reads", "unmapped reads","run name","lane","barcode"]
    #namedtuples can't have spaces
    import collections
    BamQC=collections.namedtuple('BamQC',map(lambda x: x.replace(' ','_'), fields))
    import json
    btuples=list()
    for b in bamqcs:
        f=open(b,'r')
        objects = json.load(f)
        bqc=dict((f, objects[f]) for f in fields)
        bqc2={k.replace(' ','_'): v for k,v in bqc.iteritems()}
        cur_bqc=BamQC(**bqc2)
        btuples.append(cur_bqc)
        f.close()
    return btuples

def get_files(d,wf,files):
    for k, v in d.iteritems():
        if isinstance(v, dict):
            files=get_files(v,wf,files)
        else:
            paths=set(map(lambda y : y.filePath, filter(lambda x : x.workflowName == wf, v)))
            files.add(*paths)
    return files


def read_fpr(json_file):
    # namedtuples for libs
    fields=['iusTags', 'studyTitles','sampleNames','rootSampleNames','sequencerRunNames','laneNumbers','sampleAttributes']
    import collections
    SeqLibrary=collections.namedtuple('SeqLibrary',list(fields+['is_ffpe']))
    # namedtuples for wfrs
    workflows=['FastQC','CoverageAnalysis','BamQC']
    workflow_fields=['workflowName','workflowVersion','fileSize','filePath','fileMetaType']
    WorkflowRun=collections.namedtuple('WorkflowRun',workflow_fields)
    #make collections
    from collections import defaultdict
    srli_libs=defaultdict( lambda: defaultdict(lambda: defaultdict( int )))
    srli_workflows=defaultdict( lambda: defaultdict(lambda: defaultdict( int )))
    #open that file
    import ijson
    f=open(json_file,'r')
    objects = ijson.items(f, 'item')
    for record in objects:
        sr="_".join(record['sequencerRunNames'])
        ln="_".join(record['laneNumbers'])
        it="_".join(record['iusTags'])
        if srli_libs[sr][ln][it] == 0:
            srli_libs[sr][ln]['multiplex']+=1
            isffpe=(True if is_ffpe(record) else False)
            record_fields=dict((f, record[f]) for f in fields)
            record_fields['is_ffpe']=isffpe
            cur_lib=SeqLibrary(**record_fields)
            srli_libs[sr][ln][it]=cur_lib
        if record['workflowName'] in workflows:
            if srli_workflows[sr][ln][it] == 0:
                srli_workflows[sr][ln][it]=list()
            wfr_fields=dict((f, record[f]) for f in workflow_fields)
            cur_wfr=WorkflowRun(**wfr_fields)
            srli_workflows[sr][ln][it].append(cur_wfr)
    f.close()
    return (srli_libs, srli_workflows)




def is_ffpe(record):
    if 'sampleAttributes' in record:
        if 'geo_tissue_preparation' in record['sampleAttributes']:
            for tp in record['sampleAttributes']['geo_tissue_preparation']:
                if tp == 'FFPE':
                    return True
    return False


if __name__ == "__main__":
    #assume the first file is the one we want (FPR json)
    import sys
    main(sys.argv[1])
