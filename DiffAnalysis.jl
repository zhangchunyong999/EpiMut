# EpiMut - DiffAnalysis.jl
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

using GenomicFeatures,CSV,DataFrames,Statistics,MultipleTesting,HypothesisTests,Dates,Getopt
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

function findoverlaps(query,subject)
    query_numbered= query|> enumerate .|> number_interval
    subject_numbered=subject |> enumerate .|> number_interval
    df = Vector{Tuple{Int64, Int64}}()
    for (q,r) in eachoverlap(query_numbered,subject_numbered)
        result=(
             GenomicFeatures.metadata(q).i,
             GenomicFeatures.metadata(r).i
        )
        push!(df,result)

    end
    rename!(DataFrame(df),[:queryHits,:subjectHits])
end

function generatesubject(epifile)
    a=CSV.read(epifile,DataFrame)
    kind=propertynames(a)[3]
    a1=filter(kind=>x->!isequal(x,NaN),a)
    epi=a1[!,kind]
    subject=Interval.(a1[!,:chr],a1[!,:position],a1[!,:position])
    return subject,epi
end

function classify_reads(index,match_region_cpg,epi,region_epinumber)
    epiindex = match_region_cpg[index][!,2]
    if length(epiindex)<region_epinumber
        return missing
    end
    epis=epi[epiindex]
    meanepi=mean(epis)
    return meanepi
end
function calculateregionepi(match_region_cpg,epi,region_epinumber)
    p=Vector{Union{Missing,Float64}}()
    for i in 1:match_region_cpg.ngroups
        a=classify_reads(i,match_region_cpg,epi,region_epinumber)
        push!(p,a)
    end
    p
end

function generatesumdf(overlap,query,epi,region_epinumber)
    position=unique(overlap[!,1])
    positiondf=@view DataFrame(query)[position,:]
    match_region_cpg=groupby(overlap,:1)
    regionepi=calculateregionepi(match_region_cpg,epi,region_epinumber)
    sumdf=hcat(positiondf[!,[:seqname,:first,:last]],DataFrame(pdr=regionepi))
    dropmissing!(sumdf,propertynames(sumdf)[4])
    rename!(sumdf,:seqname=>:chr,:first=>:start,:last=>:end)
    return sumdf
end
function eachfile(epifile,query,region_epinumber)
    subject,epi=generatesubject(epifile)
    overlap=findoverlaps(query,subject)
    sumdf=generatesumdf(overlap,query,epi,region_epinumber)
    return sumdf
end


function dottest(v1,v2)
    try pvalue(UnequalVarianceTTest(v1,v2))
    catch e
        missing
    end
end

function generatealldf(group,query,region_epinumber)
    dic=Dict{String,DataFrame}()
    for i in 1:group.ngroups
        group1=group[i][!,1]
        epidf=DataFrame("chr"=>String[],"start"=>Int[],"end"=>Int[])
        for i in group1
            matrix=eachfile(i,query,region_epinumber)
            epidf=outerjoin(epidf,matrix,on=[:chr,:start,:end],makeunique=true)
        end
    dic["group"*string(i)]=epidf
    end
    return innerjoin(dic["group1"],dic["group2"],on=[:chr,:start,:end],makeunique=true)
end
function generatettest(alldf,input)
    kind=unique(input[!,2])
    n1=count(x->isequal(x,kind[1]),input[!,2])
    n2=count(x->isequal(x,kind[2]),input[!,2])
    a1right=4+n1-1;a2left=a1right+1;a2right=ncol(alldf)
    m=map(eachrow(alldf)) do row
        a1 = collect(row[4:a1right])
        a2 = collect(row[a2left:a2right])
        v1=collect(skipmissing(a1));v2=collect(skipmissing(a2))
        pvalue=dottest(v1,v2)
        if isequal(pvalue,missing)
            missing,pvalue
        else
            mean(v1)-mean(v2),pvalue
        end
    end
    m
end
function generategroup(input)
    if (ncol(input)!=2) || (typeof(input[:,1])!=Vector{String}) || (length(unique(input[:,2]))!=2)
        throw(ErrorException("""The data frame you input must be 2 columns!
                        The first column must contain the correct sample directory,
                        as well as the second containing each sample's group(2 groups)!
                """))
    else
        input=sort(input,2)
        group=groupby(input,2)
    end
    group
end
function generatequery(bedfile,header,addchrname)
    if header==true
        bed=CSV.read(bedfile,DataFrame,header=1)
    else
        bed=CSV.read(bedfile,DataFrame,header=0)
    end
    if addchrname==true
        bed[!,:1]="chr".*bed[!,:1]
    end
    if (ncol(bed)<3) || (typeof(bed[:,2])!=Vector{Int}) || (typeof(bed[:,3])!=Vector{Int})
        throw(ErrorException("""The bed you input must be at least 3 columns!
                        The column must contain the chromosome name,
                        as well as the second and the third containing the position!
                """))
    else
        query=Interval.(bed[!,1],bed[!,2],bed[!,3])
    end
    sort(query)
end
function generatefinaldf(group,input,query,region_epinumber)
    alldf=generatealldf(group,query,region_epinumber)
    pvaldf=rename!(DataFrame(generatettest(alldf,input)), [1 => :delta_epimut_rate, 2 => :pval])
    h=hcat(alldf,pvaldf)
    dropmissing!(h,[:pval,:delta_epimut_rate])
    mergedf=hcat(h,MultipleTesting.adjust(h[!,:pval],BenjaminiHochberg()))
    finaldf=select(mergedf,:chr,:start,:end,:delta_epimut_rate,:pval,:x1)
    rename!(finaldf,:x1=>:qval)
    return finaldf
end
function differentialanalysis(inputfile,bedfile,header,addchrname,region_epinumber,name,outputfile)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Read the input data frame")
    input=CSV.read(inputfile,DataFrame)

    group=generategroup(input)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Read the input bed")
    query=generatequery(bedfile,header,addchrname)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Start differential analysis")
    finaldf=generatefinaldf(group,input,query,region_epinumber)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Finish differential analysis")
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Start writing the results to a csv")
    CSV.write(outputfile*'/'*name*".diffanalysis.csv",finaldf)
    println('[',Dates.format(now(), "YYYY-m-d HH:MM:SS"),']'," Finish writing the results to a csv")
end
function Argparse()
    lst = Dict{String,String}()
    for (opt, arg) in getopt(ARGS, "hi:b:H:a:r:n:o:", ["help","inputfile=","bedfile=","header=","addchrname=","region_epinumber=","name=","outputfile="])
        opt = replace(opt, "-" => "")
        arg = replace(arg, "=" => "")
        if opt == "help" || opt == "h"
            println("""
DESCRPTION      This program is for differential analysis between 2 groups using unpaired two-tailed Student's t-test.

                The input file should contain 2 columns (samplefile, group).     The bed file should contain at least 3 columns (chr, start, end).
                 Row │           samplefile                group                 Row │ chr      start     end
               ──────────────────────────────────────────────────            ─────────────────────────────────────
                   1 │ /home/zhangchunyong/group1_1.csv      1                     1 │ chr1         1      1000
                   2 │ /home/zhangchunyong/group1_2.csv      1                     2 │ chr1      1001      2000
                   3 │ /home/zhangchunyong/group2_1.csv      2                     3 │ chr2      4001      5000
                   4 │ /home/zhangchunyong/group2_2.csv      2                     4 │ chr2      5001      6000
                   5 │ /home/zhangchunyong/group1_3.csv      1                     5 │ chr10     8001      9000
                   6 │ /home/zhangchunyong/group2_3.csv      2                     6 │ chr10     9001     10000




USAGE           julia DiffAnalysis.jl [options] -i <inputfile> -b <bedfile> -n <outputname> -o <outputfile>

                    options
                    -h      --help              Show this message
                    -i      --inputfile         Input your analysis table, which contains 2 columns (samplefile, group).
                    -b      --bedfile           Input your region file (.bed).
                    -H      --header            Input whether your bed has the header. [Default: true]
                    -a      --addchrname        If the first column of your bed is "1" but not "chr1", you should input true. [Default: false]
                    -r      --region_epinumber  Input at least the number of epimut values included in each region of bed. [Default: 5]
                    -n      --name              Input output name, stored as a csv format.
                    -o      --outputfile        Input output file path.
                    -t      --threads           Input number of threads to generate differential analysis. And this parameter should follow immediately after julia,
					    	like "julia -t 20 ***.jl ...". [Default: 1]
AUTHOR
                    Contact:     Chunyong Zhang; zhangchunyong@tmu.edu.cn
""")
        elseif opt == "inputfile" || opt == "i"
            lst["inputfile"] = arg
        elseif opt == "bedfile" || opt == "b"
            lst["bedfile"] = arg
        elseif opt == "header" || opt == "H"
            lst["header"] = arg
        elseif opt == "addchrname" || opt == "a"
            lst["addchrname"] = arg
        elseif opt == "region_epinumber" || opt == "r"
            lst["region_epinumber"] = arg
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
    if isequal(c,"header") && isequal(getkey(args,c,"no"),"no")
        args[c]="true"
    elseif isequal(c,"addchrname") && isequal(getkey(args,c,"no"),"no")
        args[c]="false"
    elseif isequal(c,"region_epinumber") && isequal(getkey(args,c,"no"),"no")
        args[c]="5"
    else
        args
    end
end
args = Argparse()
if length(args) > 1
    default(args,"header");default(args,"addchrname");default(args,"region_epinumber")
    differentialanalysis(args["inputfile"], args["bedfile"], parse(Bool,args["header"]),parse(Bool,args["addchrname"]),parse(Int,args["region_epinumber"]),args["name"],args["outputfile"])
end
