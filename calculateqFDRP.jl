# EpiMut - calculateqFDRP.jl
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

using DataFrames,CSV,XAM,FASTX,GenomicFeatures,BioSequences,StatsBase,Statistics,Dates,Getopt
function f(record)
    if BAM.ispositivestrand(record)==true
        return '+'
    else
        return '-'
    end
end
function trimScigar(cigar)
	if 'S' in cigar
        cigar=replace(cigar,r"\d*S"=>s"")
    end
    return cigar
end
function trimSseq(cigar,seq)
    if 'S' in cigar
        s=count('S',cigar)
        reg=r"\d+(?=S)"
        number=[parse(Int, m.match) for m in eachmatch(reg, cigar)]
        if s==2
            seq=seq[1+number[1]:end-number[2]]
        else
            if isequal(cigar[end],'S')
                seq=seq[1:end-number[1]]
            else
                seq=seq[1+number[1]:end]
            end
        end
    end
    seq
end
function insert_at(original, i, new)
    return original[begin:i-1] * new * original[i:end]
end
function delete_at(original,left,right)
    return SubString(original,1,left-1)*SubString(original,right+1)
end

function rmI(c,reads)
    c1=replace(c,r"\d*D"=>s"")
    f=findall('I',c1)

    for i in f
        c2=c1[1:i-1]
        c2=replace(c2,r"\d*I"=>s"")
        p=parse.(Int,split(c2,'M'))
        pright=sum(p)
        pleft=sum(p[1:end-1])+1
        reads=deleteat!(reads,pleft:pright)
    end
    reads
end
function addD(c,reads)
    d1=replace(c,r"\d*I"=>s"")
    s=split.(split(d1,'D'),'M')
    len=length(s)
    right=0
    for i in 1:len-1

        p=parse.(Int,s[i])
        left=right+sum(p[1:end-1])+1
        right=right+sum(p)
        reads=insert_at(reads, left, repeat(dna"N",right-left+1))
    end
    reads
end
function indel(cigar,seq)
    #first deal with insertion,rm
    seq=rmI(cigar,seq)
    #last deal with deletion,add
    reads=addD(cigar,seq)
    return reads
end
function lastreads(cigar,reads)
    if ( 'I' in cigar) & ('D' in cigar)
        reads=indel(cigar,reads)
    elseif 'I' in cigar
        reads=rmI(cigar,reads)
    elseif 'D' in cigar
        reads=addD(cigar,reads)
    else
        reads
    end
    reads
end
function cigarreads(cigar,reads)
    reads=trimSseq(cigar,reads)
    cigar=trimScigar(cigar)
    reads=lastreads(cigar,reads)
    return cigar,reads
end
function generateref(refseq)
    ref=SubString(refseq,4,lastindex(refseq)-1)
    delposition=findfirst('_',ref)
    ref=delete_at(ref,delposition,delposition)
    return ref
end

function negmeunmeindex(reads,ref,lp,rightposition)
    ref=@view ref[lp:rightposition]
    index=findall(ExactSearchQuery(dna"CG"),LongDNA{4}(ref))
    methindex=Int[];unmethindex=Int[]
    if !isempty(index)
        for i in index
            right=i[2]-1
            onebase=reads[right]
            if isequal(onebase,DNA_G)
                push!(methindex,lp+right)
            end
            if isequal(onebase,DNA_A)
                push!(unmethindex,lp+right)
            end
        end
        if !isempty(methindex) || !isempty(unmethindex)
            return methindex,unmethindex
        end
    end
end
function posmeunmeindex(reads,ref,lp,rightposition)
    ref=@view ref[lp+1:rightposition+1]
    index=findall(ExactSearchQuery(dna"CG"),LongDNA{4}(ref))
    methindex=Int[];unmethindex=Int[]
    if !isempty(index)
        for i in index
            left=i[1]
            onebase=reads[left]
            if isequal(onebase,DNA_C)
                push!(methindex,lp+left)
            end
            if isequal(onebase,DNA_T)
                push!(unmethindex,lp+left)
            end
        end
        if !isempty(methindex) || !isempty(unmethindex)
            return methindex,unmethindex
        end
    end
end
function unmethdict(unmethCGindex)
    dic=Dict{Int64,Bool}()
    for i in unmethCGindex
        dic[i]=false
    end
    dic
end
function methdict(methCGindex)
    dic=Dict{Int64,Bool}()
    for i in methCGindex
        dic[i]=true
    end
    dic
end

function generatedic(methCGindex,dic,chromosomename,leftposition,rightposition,dmetadata)
    for i in methCGindex
        if i in keys(dic)
            push!(dic[i],Interval(chromosomename,leftposition,rightposition,'?',dmetadata))
        else
            dic[i]=[Interval(chromosomename,leftposition,rightposition,'?',dmetadata)]
        end
     end
    dic
end
function generatechrdic(bamreader,fastareader,chromosomename,final)
    dic=Dict{Int64, Vector{Interval{Dict{Int64, Bool}}}}()
    ref=extract(fastareader,chromosomename)
    for record in eachoverlap(bamreader,chromosomename,1:final)

        leftposition=BAM.position(record);rightposition=BAM.rightposition(record)
        lp=leftposition-1
        Strand=f(record)
        cigar,reads=cigarreads(BAM.cigar(record),BAM.sequence(record))
        if isequal('-',Strand)
            a=negmeunmeindex(reads,ref,lp,rightposition)
            if !isequal(a,nothing)
                methCGindex=a[1];unmethCGindex=a[2]
                dmetadata=merge(methdict(methCGindex),unmethdict(unmethCGindex))
                append!(methCGindex,unmethCGindex)
                dic=generatedic(methCGindex,dic,chromosomename,leftposition,rightposition,dmetadata)
            end
        else
            a=posmeunmeindex(reads,ref,lp,rightposition)
            if !isequal(a,nothing)
                methCGindex=a[1];unmethCGindex=a[2]
                dmetadata=merge(methdict(methCGindex),unmethdict(unmethCGindex))
                append!(methCGindex,unmethCGindex)
                dic=generatedic(methCGindex,dic,chromosomename,leftposition,rightposition,dmetadata)
            end
        end
    end
    dic
end
function number_interval(tp::Tuple)
    # Unpack Tuple.
    (i, interval) = tp

    # Setup numbered metadata.
    new_metadata = (
        i = i,
        original = GenomicFeatures.metadata(interval)
    )

    # Create new interval with numbered metadata.
    return Interval(
        seqname(interval),
        leftposition(interval),
        rightposition(interval),
        strand(interval),
        new_metadata
    )
end
function overlap_span(interval_a, interval_b)
	left = max(leftposition(interval_a), leftposition(interval_b))
	right = min(rightposition(interval_a), rightposition(interval_b))
	return length(left:right)
end

function findoverlaps(query,subject,minoverlap)
    query_numbered= query|>sort|> enumerate .|> number_interval
    subject_numbered=subject|>sort |> enumerate .|> number_interval
    itr = eachoverlap(query_numbered, subject_numbered)
    itr = Iterators.filter(itr) do (q, s)
        if overlap_span(q, s) >= minoverlap
            return true
        end
        return false
    end
    itr = Base.Generator(itr) do (q, s)
        result = (
            GenomicFeatures.metadata(q).i,
            GenomicFeatures.metadata(s).i
        )
        return result
    end
    return collect(itr)
end
function g_distance(both,sortdickey)
    a=both.-sortdickey
    return [abs(i) for i in a]
end
function generate_validinterval(i,both,sortdickey,windowsize)
    intervals=i:(i+windowsize-1)
    (maximum(intervals)>maximum(both)) && return missing
    !(sortdickey in intervals) && return missing
    (i,i+windowsize-1,count(∈(both),intervals))
end
function generate_allvalidintervals(m,both,sortdickey,windowsize)
    v=Vector{Union{Missing, Tuple{Int64, Int64, Int64}}}()
    for i in m
        gv=generate_validinterval(i,both,sortdickey,windowsize)
        push!(v,gv)
    end
    v
end
function restrict(both,sortdickey,windowsize)
    !(sortdickey in both) && return NaN
    distance=g_distance(both,sortdickey)
    both=sort(both[distance.<=windowsize])
    m=minimum(both):maximum(both)
    allvalidintervals=generate_allvalidintervals(m,both,sortdickey,windowsize)
    allvalidintervals=collect(skipmissing(allvalidintervals))
    length(allvalidintervals)==0 && return both
    maxindex=maximum(last.(allvalidintervals))
    selinterval=allvalidintervals[findfirst(x -> last(x) == maxindex,allvalidintervals)]
    return both[both.>=selinterval[1] .&& both.<=selinterval[2]]
end
function compute_discordant(index,query,subject,values,sortdickey,readpaircg,windowsize)
    q=query[index];s=subject[index]
    v1=values[q]
    length(v1)==0 && return missing
    v2=values[s]
    length(v2)==0 && return missing
    both=collect(intersect(keys(v1),keys(v2)))
    both=restrict(both,sortdickey,windowsize)
    length(both)<readpaircg && return missing
    count(x->isequal(x,NaN)||isequal(x,missing),both)>0 && return missing
    intersection1 = Union{Int64, Missing}[get(v2, k, missing) for k in both]
    intersection2=Union{Int64, Missing}[get(v1, k, missing) for k in both]
    discordant=sum(skipmissing(intersection1).!=skipmissing(intersection2))/length(intersection1)
    return discordant
end

function compute_qFDRP(chromosomename,sortdickeyvalue,minoverlap,windowsize,readpaircg,coverage)
    sortdickey=sortdickeyvalue[1];sortdicvalue=sortdickeyvalue[2]
    if length(sortdicvalue)<coverage
        return chromosomename,sortdickey,NaN,NaN,NaN
    end
    if length(sortdicvalue)>40
        sortdicvalue=sample(sortdicvalue,40,replace=false)
    end
    overlap=findoverlaps(sortdicvalue,sortdicvalue,minoverlap)
    if length(overlap)==0
        return chromosomename,sortdickey,NaN,NaN,NaN
    end
    fi=first.(overlap);la=last.(overlap)
    index=fi.<la
    values=GenomicFeatures.metadata.(sortdicvalue)
    query=fi[index];subject=la[index]
    ret=1:length(query)
    ret=[compute_discordant(i,query,subject,values,sortdickey,readpaircg,windowsize) for i in ret]
    sr=skipmissing(ret);lr=length(ret)
    isempty(sr) && return  chromosomename,sortdickey,NaN,NaN,lr
    ssr=sum(sr)
    qfdrptuple =chromosomename,sortdickey,ssr/lr,ssr,lr
    return qfdrptuple
end
function readfile(bamfile)
    if isfile(bamfile) && isfile(bamfile*".bai")
    	bamreader=open(BAM.Reader,bamfile,index=bamfile*".bai")
    else
        throw(ErrorException("There is no such a BAM/bai (the BAM index generated by samtools) file or directory"))
    end
end
function genomefile(fastafile)
    if isfile(fastafile) && isfile(fastafile*".fai")
    	fastareader=open(FASTA.Reader,fastafile,index=fastafile*".fai")
    else
        throw(ErrorException("There is no such a fasta/fai (the fasta index generated by samtools) file or directory"))
    end
end

function generateallreader(bamfile,fastafile,seqnames)
    [(readfile(bamfile),genomefile(fastafile)) for i in 1:length(seqnames)]
end

function genearteitr(bamfile,fastafile,seqnames,seqlens)
    allreaders=generateallreader(bamfile,fastafile,seqnames)
    itr=[(i,j,z) for (i,j,z) in zip(seqnames,seqlens,allreaders)]
    return itr
end

function generateqFDRPdf(chromosomename,dic,minoverlap,windowsize,readpaircg,coverage,chrfilesdict)
    sortdic=collect(sort(dic))
    qFDRP=DataFrame(compute_qFDRP.(chromosomename,sortdic,minoverlap,windowsize,readpaircg,coverage))
    rename!(qFDRP,[:chr,:position,:qFDRP,:sum_discordant_ratio,:sum_readpairs])
    chrfilesdict[chromosomename]=qFDRP
end
function generateallqFDRP(chrfilesdict,qdf)
    chrname=sort(collect(keys(chrfilesdict)))
    for chr in chrname
        append!(qdf,chrfilesdict[chr])
    end
    qdf
end
function filterchr(bamreader,rmMXY)
    name=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'M','X','Y']
    if rmMXY==true
        deleteat!(name,23:25)
    end
    allchrname=sort("chr".*string.(name))
    seqnames=bamreader.refseqnames;seqlens=bamreader.refseqlens
    index=[findfirst(isequal(i),seqnames) for i in allchrname]
    seqnames=seqnames[index];seqlens=seqlens[index]
    return seqnames,seqlens
end

function calculate_qFDRPs(bamfile::String,fastafile::String,name::String,outputfile::String,minoverlap::Int,windowsize::Int,readpaircg::Int,coverage::Int,rmMXY::Bool)
    coverage==1 && error("The coverage you input must be more than 1.")
    bamreader=readfile(bamfile)
    fastareader=genomefile(fastafile)
    seqnames,seqlens=filterchr(bamreader,rmMXY)
    itr=genearteitr(bamfile,fastafile,seqnames,seqlens)
    chrfilesdict=Dict{String,DataFrame}()
    qdf=DataFrame(chr=String[],position=Int64[],qFDRP=Float64[],sum_discordant_ratio=Float64[],sum_readpairs=Float64[])
    Threads.@threads for (chromosomename,final,bamfastareader) in itr
        println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Start calculating ",chromosomename,"'s qFDRP")
        generateqFDRPdf(chromosomename,generatechrdic(bamfastareader[1],bamfastareader[2],chromosomename,final),minoverlap,windowsize,readpaircg,coverage,chrfilesdict)
		close(bamfastareader[1]);close(bamfastareader[2])
		println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Finish calculating ",chromosomename,"'s qFDRP")
    end
    allqFDRPdf=generateallqFDRP(chrfilesdict,qdf)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Start writing all chromosomes' qFDRP results to a csv")
    CSV.write(outputfile*'/'*name*".csv",allqFDRPdf)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Finish writing all chromosomes' qFDRP results to a csv")
end
function Argparse()
    lst = Dict{String,String}()
    for (opt, arg) in getopt(ARGS, "hb:g:n:o:m:w:r:c:R:", ["help","bamfile=","genomefile=","name=","outputfile=","minoverlap=","windowsize=","readpaircg=","coverage=","rmMXY="])
        opt = replace(opt, "-" => "")
        arg = replace(arg, "=" => "")
        if opt == "help" || opt == "h"
            println("""
DESCRPTION      This program is for calculating quantative fraction of discordant read pairs (qFDRP) from the BAM.
                The BAM index (.bai) should be in the same file as the BAM. The fasta index (.fai) should be in the same file as the fasta.

USAGE           julia calculateqFDRP.jl [options] -b <in.bam> -g <in.fa> -n <outputname> -o <outputfile>

                    options
                    -h      --help          Show this message
                    -b      --bamfile       Input your BAM file, should be in the same file as the bai (generated by samtools index).
                    -g      --genomefile    Input your genome file, should be in the same file as the fai (generated by samtools faidx).
                    -m      --minoverlap    Input at least the number of overlap bases of each read pair. [Default: 35]
                    -w      --windowsize    Input the maximum distance used to restrict the concordance/discordance classification of each read pair. [Default: 50]
                    -r      --readpaircg    Input at least the number of each read pair's overlaping C loci. [Default: 2]
                    -c      --coverage      Input at least the number of reads covering a C locus.The coverage must be more than 1. [Default: 2]
                    -R      --rmMXY         When using this parameter as true, we only calculate qFDRP of autosomes without chrM, chrX, chrY.
                                            When using this parameter as false, we calculate qFDRP of all chromosomes. [Default: true]
                    -n      --name          Input output name, stored as a csv format.
                    -o      --outputfile    Input output file path.
                    -t      --threads       Input number of threads to calculate qFDRP. And this parameter should follow immediately after julia,
					    like "julia -t 20 ***.jl ...". [Default: 1]
AUTHOR
                    Contact:     Chunyong Zhang; zhangchunyong@tmu.edu.cn
""")
        elseif opt == "bamfile" || opt == "b"
            lst["bamfile"] = arg
        elseif opt == "genomefile" || opt == "g"
            lst["genomefile"] = arg
        elseif opt == "minoverlap" || opt == "m"
            lst["minoverlap"] = arg
        elseif opt == "windowsize" || opt == "w"
            lst["windowsize"] = arg
        elseif opt == "readpaircg" || opt == "r"
            lst["readpaircg"] = arg
        elseif opt == "coverage" || opt == "c"
            lst["coverage"] = arg
        elseif opt == "rmMXY" || opt == "R"
            lst["rmMXY"] = arg
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
    if isequal(c,"minoverlap") && isequal(getkey(args,c,"no"),"no")
        args[c]="35"
    elseif isequal(c,"windowsize") && isequal(getkey(args,c,"no"),"no")
        args[c]="50"
    elseif isequal(c,"readpaircg") && isequal(getkey(args,c,"no"),"no")
        args[c]="2"
    elseif isequal(c,"coverage") && isequal(getkey(args,c,"no"),"no")
        args[c]="2"
    elseif isequal(c,"rmMXY") && isequal(getkey(args,c,"no"),"no")
        args[c]="true"
    else
        args
    end
end
args = Argparse()
if length(args) > 1
    default(args,"minoverlap");default(args,"windowsize");default(args,"readpaircg");default(args,"coverage");default(args,"rmMXY")
    calculate_qFDRPs(args["bamfile"],args["genomefile"], args["name"], args["outputfile"],parse(Int,args["minoverlap"]),parse(Int,args["windowsize"]),parse(Int,args["readpaircg"]),parse(Int,args["coverage"]),parse(Bool,args["rmMXY"]))
end
