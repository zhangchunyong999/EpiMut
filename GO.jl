# EpiMut - GO.jl
#
# Copyright (C) Chunyong Zhang
# Contact: Chunyong Zhang <zhangchunyong@tmu.edu.cn>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

using RCall,Getopt
function GO(signdif,sorttss,distance,species,pvaluecutoff,qvaluecutoff,name,outputfile)
    bedtools=Sys.which("bedtools")
    run(pipeline(`$bedtools closest -a $signdif -b $sorttss -d`,
            `awk "\$NF<$distance"`,
            `awk -F'\t' '{print $(NF-1)}'`,
            `sort`,`uniq`,
            "$outputfile/$name.gene.txt"))
    if species=="hg"
        R"""
        suppressMessages(library(clusterProfiler))
        suppressMessages(library(org.Hs.eg.db))
        library(ggplot2)
        gene<-read.table(paste0($outputfile,'/',$name,'.gene.txt'))
        geneid<-mapIds(x=org.Hs.eg.db,keys=gene$V1,
                     keytype = 'SYMBOL',column = "ENTREZID")
        enrich.go.BP<-enrichGO(gene=geneid,OrgDb = org.Hs.eg.db,
                               keyType = "ENTREZID",
                               ont="BP",
                               pvalueCutoff = $pvaluecutoff,qvalueCutoff= $qvaluecutoff)
        dotplot(enrich.go.BP)
        ggsave(paste0($outputfile,'/',$name,'.enrichgo.bp.pdf'))
        write.csv(as.data.frame(enrich.go.BP),paste0($outputfile,'/',$name,'.enrichgo.bp.csv'))
        """
    else
        R"""
        suppressMessages(library(clusterProfiler))
        suppressMessages(library(org.Mm.eg.db))
        library(ggplot2)
        gene<-read.table(paste0($outputfile,'/',$name,'.gene.txt'))
        geneid<-mapIds(x=org.Mm.eg.db,keys=gene$V1,
                     keytype = 'SYMBOL',column = "ENTREZID")
        enrich.go.BP<-enrichGO(gene=geneid,OrgDb = org.Mm.eg.db,
                               keyType = "ENTREZID",
                               ont="BP",
                               pvalueCutoff = $pvaluecutoff,qvalueCutoff= $qvaluecutoff)
        dotplot(enrich.go.BP)
        ggsave(paste0($outputfile,'/',$name,'.enrichgo.bp.pdf'))
        write.csv(as.data.frame(enrich.go.BP),paste0($outputfile,'/',$name,'.enrichgo.bp.csv'))
        """
    end
end
function Argparse()
    lst = Dict{String,String}()
    for (opt, arg) in getopt(ARGS, "hi:t:d:s:p:q:n:o:", ["help","inputfile=","tssfile=","distance=","species=","pvaluecutoff=","qvaluecutoff=","name=","outputfile="])
        opt = replace(opt, "-" => "")
        arg = replace(arg, "=" => "")
        if opt == "help" || opt == "h"
            println("""
DESCRPTION      This program is for gene ontology analysis after differential analysis.
                Make sure your BEDtools was installed. We try to search your BEDtools.
                Also make sure your R has installed clusterprofiler, org.Hs.eg.db / org.Mm.eg.db, ggplot2 package.

USAGE           julia GO.jl [options] -i <inputfile> -t <tssfile> -d <distance> -s <species> -p <pvaluecutoff> -q <qvaluecutoff> -n <outputname> -o <outputfile>
​
                    options
                    -h      --help              Show this message
                    -i      --inputfile         Input your processed (signficant up or down regions) and
                                                sorted (use bedtools sort) differential analysis table (.bed).
                    -t      --tssfile           Input all genes' transcrpition start sites (TSS) table (.bed).
                    -d      --distance          Distance cutoff between significant feature and TSS. [Default: 1000]
                    -s      --speices           Species for GO analysis (hg, mm). [Default: hg]
                    -p      --pvaluecutoff      Adjusted pvalue cutoff on enrichment tests to report. [Default: 0.05]
                    -q      --qvaluecutoff      Qvalue cutoff on enrichment tests to report as significant. [Default: 0.2]
                    -n      --name              Input output name, stored as a pdf and csv format.
                    -o      --outputfile        Input output file path.

AUTHOR
                    Contact:     Chunyong Zhang; zhangchunyong@tmu.edu.cn
""")
        elseif opt == "inputfile" || opt == "i"
            lst["inputfile"] = arg
        elseif opt == "tssfile" || opt == "t"
            lst["tssfile"] = arg
        elseif opt == "distance" || opt == "d"
            lst["distance"] = arg
        elseif opt == "species" || opt == "s"
            lst["species"] = arg
        elseif opt == "pvaluecutoff" || opt == "p"
            lst["pvaluecutoff"] = arg
        elseif opt == "qvaluecutoff" || opt == "q"
            lst["qvaluecutoff"] = arg
        elseif opt == "name" || opt == "n"
            lst["name"] = arg
        elseif opt == "outputfile" || opt == "o"
            lst["outputfile"] = arg
        else
            println("Please check your parameter!")
        end
    end
    lst
end
function default(args,c)
    if isequal(c,"distance") && isequal(getkey(args,c,"no"),"no")
        args[c]="1000"
    elseif isequal(c,"species") && isequal(getkey(args,c,"no"),"no")
        args[c]="hg"
    elseif isequal(c,"pvaluecutoff") && isequal(getkey(args,c,"no"),"no")
        args[c]="0.05"
    elseif isequal(c,"qvaluecutoff") && isequal(getkey(args,c,"no"),"no")
        args[c]="0.2"
    else
        args
    end
end
args = Argparse()
if length(args) > 1
    default(args,"distance");default(args,"species");default(args,"pvaluecutoff");default(args,"qvaluecutoff")
    GO(args["inputfile"],args["tssfile"],parse(Int,args["distance"]),args["species"],parse(Float64,args["pvaluecutoff"]),parse(Float64,args["qvaluecutoff"]),args["name"],args["outputfile"])
end
