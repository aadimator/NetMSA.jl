module NetMSA

using StatsBase

export createPeerMatrix, matrixalignment

mutable struct Position
    row::Int64
    indexes::Vector{Int64}
end

mutable struct Particle
    value::Char
    updated::Int64
    pos::Position
    best::Position
    bestvalue::Float64

    function Particle(value::Char, pos::Position)
        return new(value, 0, pos, pos, 0.0)
    end
end

"""
    createPeerMatrix(inputStrings::Vector{String})::Matrix{Union{Missing,Char}}

Create and return a Peer matrix, containing charachters as elements, where each input sequence,
provided in the inputStrings, is represented as a column. Missing values are represented in the
matrix by the `missing` keyword.

# Examples
```jldoctest
julia> createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'
```
"""
function createPeerMatrix(inputStrings::Vector{String})::Matrix{Union{Missing,Char}}
		col_size = size(inputStrings, 1);
		row_size = max([length(s) for s in inputStrings]...);
		peer_matrix = Matrix{Union{Missing,Char}}(missing, row_size, col_size);
		for (i, s) in enumerate(inputStrings)
			peer_matrix[1:length(s), i] = collect(s);
		end
		return peer_matrix;
end


"""
    getposition(value::T, rowindex::Int64, matrix::Matrix{T}) where (T)

Return the Position (rowindex, [colindex1, colindex2, ...]) of the Particle
represented by `value`, at the `rowindex` in the `matrix`.

# Examples
```jldoctest
julia> M = createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

 juila> NetMSA.getposition('b', 2, M)
 NetMSA.Position(2, [1, 3, 4])
```
"""
function getposition(value::T, rowindex::Int64, matrix::Matrix{T}) where (T)
  indexes = findall(i -> i == value, skipmissing(matrix[rowindex, :]))
  return Position(rowindex, indexes);
end

function mostfrequent(row)
  counts = countmap(row);
  delete!(counts, '-');
  max = isempty(counts) ? 0 : findmax(counts);
  return max;
end

function aligned(row)::Bool
  row = Set(row)
  return (length(row) == 1 && !(missing in row)) || (length(row) == 2 && ('-' in row || missing in row))
end

function full(row)::Bool
  return length(Set(row)) == 1 && !(missing in Set(row))
end

function weight(row; w1=0.25, w2=0.5, w3=1.0)
  if length(Set(skipmissing(row))) == 0
    return 0;
  end
  if full(row)
    return w3;
  end

  c = length(row);
  if c == 0
    return 0;
  end
  max = mostfrequent(row)[1];

  if aligned(row)
    return w2 * max / c;
  else
    x = max <= 1 ? 0 : max;
    return w1 * x / c;
  end
end

function objective(M, rowindex::Int; endindex::Int=0)
  weights = sum(weight.(eachrow(M[rowindex:end, :])))
  C = mostfrequent(M[rowindex, :])[1];
  A = sum(aligned.(eachrow(M))[rowindex:end])

  endindex = endindex == 0 ? size(M)[1] : endindex;
  if endindex > size(M)[1]
    throw(ArgumentError("endind exceeds the matrix size"));
  end

  counts = countmap(M[rowindex:endindex, :]);
  Gaps = get(counts, '-', 0);

  return weights * (A * C) / (1 + Gaps)
end

function createswarm(rowindex::Int64, row)
  unique = Set(skipmissing(row))
  swarm = Vector{Particle}(undef, length(unique))
  for (i, c) in enumerate(unique)
    swarm[i] = Particle(c, getposition(rowindex, row, c))
  end
  return swarm
end

function criteria3(p::Particle, M, newindex)
#   display(newindex)
#   display(M[newindex, :])
#   display(length(p.pos.indexes) != length(getposition(newindex, M[newindex, :], p.value).indexes))
  return length(p.pos.indexes) != length(getposition(newindex, M[newindex, :], p.value).indexes)
end

function criteria2(p::Particle)
  return p.updated > 6;
end

function stopcriteria(p::Particle, M, t)
  c3 = criteria3(p, M, t);
  c2 = criteria2(p);
  if c3
    display("Terminating cause of criteria 3")
  elseif c2
    display("Terminating cause of criteria 2")
  end
  return c3 || c2;
end

function remove_missing_rows(M)
  return M[[length(Set(skipmissing(r))) != 0 for r in eachrow(M)], :]
end

function flydown(p, M; stride=1)
  notpcols = setdiff(collect(1:size(M, 2)), p.pos.indexes)
  colsize = size(M, 2)
  pos = p.pos
  newrows = fill('-', (stride, colsize))
  M = vcat(M[1:pos.row - 1, :], reshape(newrows, stride, colsize), M[pos.row:end, :])
  display(M)
  for i in collect(pos.row + stride:size(M, 1))
    M[i - stride,notpcols] = M[i, notpcols]
    M[i, notpcols] .= missing
  end
  display(M)
  M = remove_missing_rows(M)
  return M
end

function rowalignment(r, M)
  row = M[r, :];
#   display(row);
#   println(aligned(row))
  if aligned(row)
#     println("aligned");
    return nothing;
  end

  swarm = createswarm(r, row);

  gₒ = g = swarm[1];
  gₒvalue = gvalue = objective(M, r, endindex=r);
  cols = size(row, 1)

  for p in swarm

    t = r;
    N = copy(M);

    p.bestvalue = objective(M, r, endindex=t)
#     display("Aligning $p");
#     display(N)
#     display(p)

    missingp = setdiff(collect(1:size(N, 2)), p.pos.indexes)
    maxlen = maximum([length(collect(skipmissing(col))) for col in eachcol(N[:, missingp])])
    criteria1 = maxlen;

#     display(criteria1)

    while stopcriteria(p, N, t) != true && t < criteria1
#       display(stopcriteria(p, N, t) != true)
      t += 1;
      p.updated += 1;

      N = flydown(p, N);
#       display(N)
      display(p)
      score = objective(N, r);
      display(score);
#       display(p.bestvalue);
      if score > p.bestvalue
        p.bestvalue = score;
        p.updated = 0;
      end

      if score > gvalue
        gvalue = score;
        g = deepcopy(p);
        g.best = getposition(t, N[t, :], p.value);
        g.bestvalue = score;
      end


#       display(p)

    end
  end

  if gvalue == gₒvalue
    return nothing;
  end
  return g;
end

function matrixalignment(M)
  for (index, row) in enumerate(eachrow(M))
  #   println("$index: $row")
    g = rowalignment(index, M)
    if !isnothing(g)
      M = flydown(g, M, stride=g.best.row - g.pos.row)
    end
  end
  replace!(M, missing => '-')
end

end
