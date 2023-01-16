# EpiMut - calculateMCR.jl
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

using DataFrames,CSV,XAM,FASTX,GenomicFeatures,BioSequences,Statistics,Dates,Getopt
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
function negmeunmeindex(reads,ref,lp)
    ref=reverse_complement(LongDNA{4}(ref))
    index=findall(ExactSearchQuery(dna"CG"),ref)
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
function posmeunmeindex(reads,ref,lp)
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
function conver(methCGindex,unmethCGindex)
    lm=length(methCGindex);lum=length(unmethCGindex)
    if (lm==0) | (lum==0)
        return false,[lm,lum]
    else
        return false,true,[lm,lum]
    end
end
function generatedic1(methCGindex,dic)
    for i in methCGindex
        if i in keys(dic)
            push!(dic[i],false)
        else
            dic[i]=[false]
        end
     end
    dic
end
function generatedic2(methCGindex,dic,ismixed)
    methnum=ismixed[3][1]
    for i in 1:methnum
        m=methCGindex[i]
        if m in keys(dic)
            push!(dic[m],false)
        else
            dic[m]=[false]
        end
     end
    for j in methnum+1:length(methCGindex)
        n=methCGindex[j]
        if n in keys(dic)
            push!(dic[n],true)
        else
            dic[n]=[true]
        end
     end
    dic
end
function generatefinaldic(methCGindex,unmethCGindex,dic)
    ismixed=conver(methCGindex,unmethCGindex)
    if length(ismixed)==2
        append!(methCGindex,unmethCGindex)
        dic=generatedic1(methCGindex,dic)
    else
        append!(methCGindex,unmethCGindex)
        dic=generatedic2(methCGindex,dic,ismixed)
    end
    dic
end
function generatechrdic(bamreader,fastareader,chromosomename,final)
    dic=Dict{Int64,Vector{Bool}}()
    ref=extract(fastareader,chromosomename)
    for record in eachoverlap(bamreader,chromosomename,1:final)
        leftposition=BAM.position(record)
        lp=leftposition-1
        rightposition=BAM.rightposition(record)
        Strand=f(record)
        cigar,reads=cigarreads(BAM.cigar(record),BAM.sequence(record))
        if isequal('-',Strand)
            a=negmeunmeindex(reads,ref,lp,rightposition)
            if !isequal(a,nothing)
                methCGindex=a[1];unmethCGindex=a[2]
                dic=generatefinaldic(methCGindex,unmethCGindex,dic)
            end
        else
            a=posmeunmeindex(reads,ref,lp,rightposition)
            if !isequal(a,nothing)
                methCGindex=a[1];unmethCGindex=a[2]
                dic=generatefinaldic(methCGindex,unmethCGindex,dic)
            end
        end
    end
    dic
end

function compute_MCR(index,coverage)#weighted methylation concurrence ratio
    mconcurrence=count(x->isequal(x,true),index)
    sum=length(index)
    if sum<coverage
        return NaN,mconcurrence,sum
    else
        return mean(index),mconcurrence,sum
    end
end
function calculateMCR(sortdicvalues,coverage)#weighted methylation concurrence ratio
    MCR=Vector{Tuple{Float64,Int,Int}}()
    for index in sortdicvalues
        push!(MCR,compute_MCR(index,coverage))
    end
    MCR
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
function generateMCRdf(chromosomename,dic,coverage,chrfilesdict)
    sortdic=sort(dic)
    sortdickeys=collect(keys(sortdic))
    sortdicvalues=collect(values(sortdic))
    MCRs=calculateMCR(sortdicvalues,coverage)
    MCR=hcat(DataFrame(x=repeat([chromosomename],length(MCRs)),y=sortdickeys),DataFrame(MCRs))
    rename!(MCR,[:chr,:position,:MCR,:concurrence,:sum])
    chrfilesdict[chromosomename]=MCR
end
function generateallMCR(chrfilesdict,cdf)
    chrname=sort(collect(keys(chrfilesdict)))
    for chr in chrname
        append!(cdf,chrfilesdict[chr])
    end
    cdf
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

function calculate_MCRs(bamfile::String,fastafile::String,name::String,outputfile::String,coverage::Int,rmMXY::Bool)
    bamreader=readfile(bamfile)
    fastareader=genomefile(fastafile)
    seqnames,seqlens=filterchr(bamreader,rmMXY)
    itr=genearteitr(bamfile,fastafile,seqnames,seqlens)
    chrfilesdict=Dict{String,DataFrame}()
    cdf=DataFrame(chr=String[],position=Int64[],MCR=Float64[],concurrence=Int[],sum=Int[])
    Threads.@threads for (chromosomename,final,bamfastareader) in itr
        println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Start calculating ",chromosomename,"'s MCR")
		generateMCRdf(chromosomename,generatechrdic(bamfastareader[1],bamfastareader[2],chromosomename,final),coverage,chrfilesdict)
		close(bamfastareader[1]);close(bamfastareader[2])
		println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Finish calculating ",chromosomename,"'s MCR")
end
    allMCRdf=generateallMCR(chrfilesdict,cdf)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Start writing all chromosomes' MCR results to a csv")
    CSV.write(outputfile*'/'*name*".csv",allMCRdf)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Finish writing all chromosomes' MCR results to a csv")
end
function Argparse()
    lst = Dict{String,String}()
    for (opt, arg) in getopt(ARGS, "hb:g:n:o:c:R:", ["help","bamfile=","genomefile=","name=","outputfile=","coverage=","rmMXY="])
        opt = replace(opt, "-" => "")
        arg = replace(arg, "=" => "")
        if opt == "help" || opt == "h"
            println("""
DESCRPTION      This program is for calculating methylation concurrence ratio (MCR) from the BAM.
                The BAM index (.bai) should be in the same file as the BAM. The fasta index (.fai) should be in the same file as the fasta.

USAGE           julia calculateMCR.jl [options] -b <in.bam> -g <in.fa> -n <outputname> -o <outputfile>

                    options
                    -h      --help          Show this message
                    -b      --bamfile       Input your BAM file, should be in the same file as the bai (generated by samtools index).
                    -g      --genomefile    Input your genome file, should be in the same file as the fai (generated by samtools faidx).
                    -c      --coverage      Input at least the number of reads covering a C locus. [Default: 1]
                    -R      --rmMXY         When using this parameter as true, we only calculate MCR of autosomes without chrM, chrX, chrY.
                                            When using this parameter as false, we calculate MCR of all chromosomes. [Default: true]
                    -n      --name          Input output name, stored as a csv format.
                    -o      --outputfile    Input output file path.
                    -t      --threads       Input number of threads to calculate MCR. And this parameter should follow immediately after julia,
					    like "julia -t 20 ***.jl ...". [Default: 1]
AUTHOR
                    Contact:     Chunyong Zhang; zhangchunyong@tmu.edu.cn
""")
        elseif opt == "bamfile" || opt == "b"
            lst["bamfile"] = arg
        elseif opt == "genomefile" || opt == "g"
            lst["genomefile"] = arg
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
    if isequal(c,"coverage") && isequal(getkey(args,c,"no"),"no")
        args[c]="1"
    elseif isequal(c,"rmMXY") && isequal(getkey(args,c,"no"),"no")
        args[c]="true"
    else
        args
    end
end
args = Argparse()
if length(args) > 1
    default(args,"coverage");default(args,"rmMXY")
    calculate_MCRs(args["bamfile"],args["genomefile"], args["name"], args["outputfile"],parse(Int,args["coverage"]),parse(Bool,args["rmMXY"]))
end
