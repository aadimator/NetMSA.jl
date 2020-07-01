var documenterSearchIndex = {"docs":
[{"location":"#","page":"Home","title":"Home","text":"CurrentModule = NetMSA","category":"page"},{"location":"#NetMSA-1","page":"Home","title":"NetMSA","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [NetMSA]","category":"page"},{"location":"#NetMSA.NetMSA","page":"Home","title":"NetMSA.NetMSA","text":"This module provides an implementation of NetMSA algorithm in Julia, which can be used for multiple sequence alignment.\n\n\n\n\n\n","category":"module"},{"location":"#NetMSA.createPeerMatrix-Tuple{Array{String,1}}","page":"Home","title":"NetMSA.createPeerMatrix","text":"createPeerMatrix(inputStrings::Vector{String})::Matrix{Union{Missing,Char}}\n\nCreate and return a Peer matrix, containing charachters as elements, where each input sequence, provided in the inputStrings, is represented as a column. Missing values are represented in the matrix by the missing keyword.\n\nExamples\n\njulia> NetMSA.createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.Particle","page":"Home","title":"NetMSA.Particle","text":"A particle that is used for creating swarms.\n\nFields\n\nvalue::Char : Value of the particle, e.g. 'b' or 'c'\nupdated::Int64 : Number of turns till last updated\npos::Position : The original position of the particle\nbest::Position : The best local position of the particle\nbestvalue::Float64 : Best local score\n\n\n\n\n\n","category":"type"},{"location":"#NetMSA.Position","page":"Home","title":"NetMSA.Position","text":"Store the position of a given particle. Position x_s_i (r) of the particle p_s_i is defined by using the row r that contains the symbol s_i as well as locations of the symbol s_i in the different columns (indexes of the columns that contain s_i ) in the row r.\n\n\n\n\n\n","category":"type"},{"location":"#NetMSA.aligned-Tuple{Any}","page":"Home","title":"NetMSA.aligned","text":"aligned(row)::Bool\n\nReturn whether a row is aligned or not.\n\nA row is aligned if it only contains different occurrences of the same symbol.\n\nExamples\n\njulia> M = NetMSA.createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.aligned(M[1, :])\n true\n\n juila> NetMSA.aligned(M[2, :])\n false\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.full-Tuple{Any}","page":"Home","title":"NetMSA.full","text":"aligned(row)::Bool\n\nReturn whether a row is full or not.\n\nAn aligned row r is called full if no gaps (—) are added in the row r . That is, the number of occurrences of the symbol in the row is equal to the number of columns in the matrix.\n\nExamples\n\njulia> M = NetMSA.createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.full(M[1, :])\n true\n\n juila> NetMSA.full(M[2, :])\n false\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.getposition-Tuple{Any,Any,Any}","page":"Home","title":"NetMSA.getposition","text":"getposition(value, rowindex, matrix)\n\nReturn the Position (rowindex, [colindex1, colindex2, ...]) of the Particle represented by value, at the rowindex in the matrix.\n\nExamples\n\njulia> M = NetMSA.createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.getposition('b', 2, M)\n NetMSA.Position(2, [1, 3, 4])\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.mostfrequent-Tuple{Any}","page":"Home","title":"NetMSA.mostfrequent","text":"mostfrequent(row)\n\nReturn a tuple containing the most frequent element occuring in the row, along with its frequency.\n\nExamples\n\njulia> M = NetMSA.createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.mostfrequent(M[2, :])\n (3, 'b')\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.objective-Tuple{Any,Any}","page":"Home","title":"NetMSA.objective","text":"objective(M, rowindex; endindex=0)\n\nReturn the objective score of the row, calculated as follows:\n\nf(x_s(r)) = fracA(r) times C(r)1 + Gaps(r) times  sum_j=r^k w(j)\n\nwhere A(r) is the number of aligned rows in M from r to the last row, C(r) is the maximum number of matched charachters in the current row, Gaps(r) is the number of gaps added to the matrix M from row r to the last row, and w(r) is the weight of the row r.\n\nendindex is used to reduce the search area for Gaps, and if it is not provided, it would default to size(M)[1].\n\nExamples\n\njulia> M = NetMSA.createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.objective(M, 2)\n 2.625\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.weight-Tuple{Any}","page":"Home","title":"NetMSA.weight","text":"weight(row; w1=0.25, w2=0.5, w3=1.0)\n\nReturn the weight of the row, calculated as:\n\nw(r) = begincases\n  w_1  times  fracxc  text if r is not aligned   \n  w_2  times  fracn_sc  text if r is aligned   \n  w_3  text if r is full   \nendcases\n\nwhere n_s is the number of occurrences of the symbol s in the aligned row r, and c is the total number of columns in the row. The value of x is equal to zero if every symbol in the row r occurred at most once, otherwise x is equal to the max number of occurrences (matches) of some symbol in r.\n\nExamples\n\njulia> M = NetMSA.createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.weight(M[1, :])\n 1.0\n\n juila> NetMSA.weight(M[2, :])\n 0.1875\n\n\n\n\n\n","category":"method"}]
}
