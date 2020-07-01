var documenterSearchIndex = {"docs":
[{"location":"#","page":"Home","title":"Home","text":"CurrentModule = NetMSA","category":"page"},{"location":"#NetMSA-1","page":"Home","title":"NetMSA","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [NetMSA]","category":"page"},{"location":"#NetMSA.createPeerMatrix-Tuple{Array{String,1}}","page":"Home","title":"NetMSA.createPeerMatrix","text":"createPeerMatrix(inputStrings::Vector{String})::Matrix{Union{Missing,Char}}\n\nCreate and return a Peer matrix, containing charachters as elements, where each input sequence, provided in the inputStrings, is represented as a column. Missing values are represented in the matrix by the missing keyword.\n\nExamples\n\njulia> createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.aligned-Tuple{Any}","page":"Home","title":"NetMSA.aligned","text":"aligned(row)::Bool\n\nReturn whether a row is aligned or not.\n\nA row is aligned if it only contains different occurrences of the same symbol.\n\nExamples\n\njulia> M = createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.aligned(M[1, :])\n true\n\n juila> NetMSA.aligned(M[2, :])\n false\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.full-Tuple{Any}","page":"Home","title":"NetMSA.full","text":"aligned(row)::Bool\n\nReturn whether a row is full or not.\n\nAn aligned row r is called full if no gaps (—) are added in the row r . That is, the number of occurrences of the symbol in the row is equal to the number of columns in the matrix.\n\nExamples\n\njulia> M = createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.full(M[1, :])\n true\n\n juila> NetMSA.full(M[2, :])\n false\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.getposition-Union{Tuple{T}, Tuple{T,Int64,Array{T,2}}} where T","page":"Home","title":"NetMSA.getposition","text":"getposition(value::T, rowindex::Int64, matrix::Matrix{T}) where (T)\n\nReturn the Position (rowindex, [colindex1, colindex2, ...]) of the Particle represented by value, at the rowindex in the matrix.\n\nExamples\n\njulia> M = createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.getposition('b', 2, M)\n NetMSA.Position(2, [1, 3, 4])\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.mostfrequent-Tuple{Any}","page":"Home","title":"NetMSA.mostfrequent","text":"mostfrequent(row)\n\nReturn a tuple containing the most frequent element occuring in the row, along with its frequency.\n\nExamples\n\njulia> M = createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.mostfrequent(M[2, :])\n (3, 'b')\n\n\n\n\n\n","category":"method"},{"location":"#NetMSA.weight-Tuple{Any}","page":"Home","title":"NetMSA.weight","text":"weight(row; w1=0.25, w2=0.5, w3=1.0)\n\nReturn the weight of the row, calculated as:\n\nw(r) =  eginequation*\negincases\n  w_1  \times  \fracxc  \text if r is not aligned   \n  w_2  \times  \fracn_sc  \text if r is aligned   \n  w_3  \text if r is full   \nndcases\nndequation*\n\nwhere n_s is the number of occurrences of the symbol s in the aligned row r, and c is the total number of columns in the row. The value of x is equal to zero if every symbol in the row r occurred at most once, otherwise x is equal to the max number of occurrences (matches) of some symbol in r.\n\nExamples\n\njulia> M = createPeerMatrix([\"abcbcdem\", \"acbcfg\", \"abchimn\", \"abcbcjkm\"])\n8×4 Array{Union{Missing, Char},2}:\n 'a'  'a'      'a'      'a'\n 'b'  'c'      'b'      'b'\n 'c'  'b'      'c'      'c'\n 'b'  'c'      'h'      'b'\n 'c'  'f'      'i'      'c'\n 'd'  'g'      'm'      'j'\n 'e'  missing  'n'      'k'\n 'm'  missing  missing  'm'\n\n juila> NetMSA.weight(M[1, :])\n 1.0\n\n juila> NetMSA.weight(M[2, :])\n 0.1875\n\n\n\n\n\n","category":"method"}]
}
