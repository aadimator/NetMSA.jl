"""
This module provides an implementation of NetMSA algorithm in Julia, which
can be used for multiple sequence alignment.
"""
module NetMSA

using StatsBase

export createPeerMatrix, matrixalignment

@doc raw"""
Store the position of a given particle. Position ``x_{s_{i}}(r)`` of the particle
``p_{s_{i}}`` is defined by using the row ``r`` that contains the symbol ``s_i`` as
well as locations of the symbol ``s_i`` in the different columns (indexes of the
columns that contain ``s_i`` ) in the row ``r``.
"""
mutable struct Position
    row::Int64
    indexes::Vector{Int64}
end

"""
A particle that is used for creating swarms.

# Fields
- value::Char : Value of the particle, e.g. 'b' or 'c'
- updated::Int64 : Number of turns till last updated
- pos::Position : The original position of the particle
- best::Position : The best local position of the particle
- bestvalue::Float64 : Best local score
"""
mutable struct Particle
    value::Char
    updated::Int64
    pos::Position
    best::Position
    bestvalue::Float64

    """
        Particle(value::Char, pos::Position)

    Initializes the Particle, where the best position is intialized
    to be the current position.
    """
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
julia> NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
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
    getposition(value, rowindex, matrix)

Return the Position (rowindex, [colindex1, colindex2, ...]) of the Particle
represented by `value`, at the `rowindex` in the `matrix`.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
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
function getposition(value, rowindex, matrix)
  indexes = findall(i -> i == value, skipmissing(matrix[rowindex, :]))
  return Position(rowindex, indexes);
end

"""
    mostfrequent(row)

Return a tuple containing the most frequent element occuring in the `row`, along
with its frequency.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

 juila> NetMSA.mostfrequent(M[2, :])
 (3, 'b')
```
"""
function mostfrequent(row)
  counts = countmap(row);
  delete!(counts, '-');
  max = isempty(counts) ? 0 : findmax(counts);
  return max;
end

"""
    aligned(row)::Bool

Return whether a row is aligned or not.

A row is *aligned* if it only contains different occurrences of the same symbol.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

 juila> NetMSA.aligned(M[1, :])
 true

 juila> NetMSA.aligned(M[2, :])
 false
```
"""
function aligned(row)::Bool
  row = Set(row)
  return (length(row) == 1 && !(missing in row)) || (length(row) == 2 && ('-' in row || missing in row))
end

"""
    aligned(row)::Bool

Return whether a row is full or not.

An aligned row r is called full if no gaps (—) are added in the row r .
That is, the number of occurrences of the symbol in the row is equal to
the number of columns in the matrix.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

 juila> NetMSA.full(M[1, :])
 true

 juila> NetMSA.full(M[2, :])
 false
```
"""
function full(row)::Bool
  return length(Set(row)) == 1 && !(missing in Set(row))
end

@doc raw"""
    weight(row; w1=0.25, w2=0.5, w3=1.0)

Return the weight of the row, calculated as:

```math
w(r) = \begin{cases}
  w_1  \times  \frac{x}{c}; & \text{ if r is not aligned}   \\
  w_2  \times  \frac{n_s}{c}; & \text{ if r is aligned}   \\
  w_3; & \text{ if r is full}   \\
\end{cases}
```
where ``n_s`` is the number of occurrences of the symbol ``s`` in the aligned row ``r``,
and ``c`` is the total number of columns in the row. The value of ``x`` is equal to zero
if every symbol in the row ``r`` occurred at most once, otherwise ``x`` is equal to the
max number of occurrences (matches) of some symbol in ``r``.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

 juila> NetMSA.weight(M[1, :])
 1.0

 juila> NetMSA.weight(M[2, :])
 0.1875
```
"""
function weight(row; w1=0.25, w2=0.5, w3=1.0)
  if full(row)
    return w3;
    end

  c = length(row);
  max = mostfrequent(row)[1];

  if aligned(row)
    return w2 * max / c;
  else
    x = max <= 1 ? 0 : max;
    return w1 * x / c;
  end
end

@doc raw"""
    objective(M, rowindex; endindex=0)

Return the objective score of the row, calculated as follows:

```math
f(x_s(r)) = \frac{A(r) \times C(r)}{1 + Gaps(r)} \times  \sum_{j=r}^{k} w(j)
```
where ``A(r)`` is the number of aligned rows in M from ``r`` to the last row,
``C(r)`` is the maximum number of matched charachters in the current row,
``Gaps(r)`` is the number of gaps added to the matrix M from row ``r`` to the
last row, and ``w(r)`` is the weight of the row ``r``.


`endindex` is used to reduce the search area for Gaps, and if it is not provided,
it would default to `size(M)[1]`.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

 juila> NetMSA.objective(M, 2)
 2.625
```
"""
function objective(M, rowindex; endindex=0)
  weights = sum(weight.(eachrow(M[rowindex:end, :])))
  C = mostfrequent(M[rowindex, :])[1];
  A = sum(aligned.(eachrow(M))[rowindex:end])

  endindex = endindex == 0 ? size(M)[1] : endindex;
  if endindex > size(M)[1]
    throw(ArgumentError("endindex exceeds the matrix size"));
    end

  counts = countmap(M[rowindex:endindex, :]);
  Gaps = get(counts, '-', 0);

  return weights * (A * C) / (1 + Gaps)
end

"""
    createswarm(rowindex::Int64, M)

Create a swarm containing unique Particles in the current row.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

juila> NetMSA.createswarm(2, M)
2-element Array{NetMSA.Particle,1}:
 NetMSA.Particle('c', 0, NetMSA.Position(2, [2]), NetMSA.Position(2, [2]), 0.0)
 NetMSA.Particle('b', 0, NetMSA.Position(2, [1, 3, 4]), NetMSA.Position(2, [1, 3, 4]), 0.0)
```
"""
function createswarm(rowindex::Int64, M)
  unique = Set(skipmissing(M[rowindex, :]))
  swarm = Vector{Particle}(undef, length(unique))
  for (i, c) in enumerate(unique)
    swarm[i] = Particle(c, getposition(c, rowindex, M))
  end
  return swarm
end

function criteria3(p::Particle, newindex::Int, M)
  return length(p.pos.indexes) != length(getposition(p.value, newindex, M).indexes)
end

function criteria2(p::Particle; threshold::Int=5)
return p.updated > threshold;
end

"""
    stopcriteria(p::Particle, newindex, M; threshold::Int=5, debug=false)

Check whether a stopping criteria has been met. Two stopping criteria
are checked in this function:
  - Criteria 2: If a particle hasn't updated its best score in the last `threshold` turns.
  - Criteria 3: If a particle moves to a new row which already contains the same symbol as that of the particle.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

juila> p = NetMSA.Particle('b', NetMSA.getposition('b', 2, M));
NetMSA.Particle('b', 0, Main.NetMSA.Position(2, [1, 3, 4]), Main.NetMSA.Position(2, [1, 3, 4]), 0.0)

julia> NetMSA.stopcriteria(p, 3, M; debug=true)
"Terminating because of criteria 3"
true
```
"""
function stopcriteria(p::Particle, newindex, M; threshold::Int=5, debug=false)
  c3 = criteria3(p, newindex, M);
  c2 = criteria2(p; threshold=threshold);
  if debug && c3
    display("Terminating because of criteria 3")
  elseif debug && c2
    display("Terminating because of criteria 2")
  end
  return c3 || c2;
end

function remove_missing_rows(M)
  return M[[length(Set(skipmissing(r))) != 0 for r in eachrow(M)], :]
end

"""
    flydown(p, M; stride=1)

Fly down the given particle by `stride`.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

juila> p = NetMSA.Particle('b', NetMSA.getposition('b', 2, M));
NetMSA.Particle('b', 0, Main.NetMSA.Position(2, [1, 3, 4]), Main.NetMSA.Position(2, [1, 3, 4]), 0.0)

julia> NetMSA.flydown(p, M)
9×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 '-'  'c'      '-'      '-'
 'b'  'b'      'b'      'b'
 'c'  'c'      'c'      'c'
 'b'  'f'      'h'      'b'
 'c'  'g'      'i'      'c'
 'd'  missing  'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

julia> NetMSA.flydown(p, M; stride=3)
11×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 '-'  'c'      '-'      '-'
 '-'  'b'      '-'      '-'
 '-'  'c'      '-'      '-'
 'b'  'f'      'b'      'b'
 'c'  'g'      'c'      'c'
 'b'  missing  'h'      'b'
 'c'  missing  'i'      'c'
 'd'  missing  'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'
```
"""
function flydown(p, M; stride=1)
  notpcols = setdiff(collect(1:size(M, 2)), p.pos.indexes)
  colsize = size(M, 2)
  pos = p.pos
  newrows = fill('-', (stride, colsize))
  M = vcat(M[1:pos.row - 1, :], reshape(newrows, stride, colsize), M[pos.row:end, :])
  for i in collect(pos.row + stride:size(M, 1))
    M[i - stride,notpcols] = M[i, notpcols]
    M[i, notpcols] .= missing
  end
  M = remove_missing_rows(M)
end

"""
    rowalignment(rowindex, M)

Return Particle with best position that aligns (maximizes the objective score)
the given row in the matrix.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

juila> p = NetMSA.Particle('b', NetMSA.getposition('b', 2, M));
NetMSA.Particle('b', 0, Main.NetMSA.Position(2, [1, 3, 4]), Main.NetMSA.Position(2, [1, 3, 4]), 0.0)

julia> NetMSA.rowalignment(2, M)
NetMSA.Particle('c', 0, Main.NetMSA.Position(2, [2]), Main.NetMSA.Position(3, [1]), 9.0)
```
"""
function rowalignment(rowindex, M)
  row = M[rowindex, :];
  if aligned(row)
    return nothing;
  end

  swarm = createswarm(rowindex, M);

  gₒ = g = swarm[1];
  gₒvalue = gvalue = objective(M, rowindex, endindex=rowindex);
  cols = size(row, 1)

  for p in swarm

    t = rowindex;
    N = copy(M);

    p.bestvalue = objective(M, rowindex, endindex=t)

    missingp = setdiff(collect(1:size(N, 2)), p.pos.indexes)
    maxlen = maximum([length(collect(skipmissing(col))) for col in eachcol(N[:, missingp])])
    criteria1 = maxlen;


    while stopcriteria(p, t, N) != true && t < criteria1
      t += 1;
      p.updated += 1;

      N = flydown(p, N);
      score = objective(N, rowindex);
      if score > p.bestvalue
        p.bestvalue = score;
        p.updated = 0;
      end

      if score > gvalue
        gvalue = score;
        g = deepcopy(p);
        g.best = getposition(p.value, t, N);
        g.bestvalue = score;
    end

    end
  end

  if gvalue == gₒvalue
    return nothing;
  end

  return g;
end


"""
    matrixalignment(M)

Align the matrix using NetMSA algorithm.

# Examples
```jldoctest
julia> M = NetMSA.createPeerMatrix(["abcbcdem", "acbcfg", "abchimn", "abcbcjkm"])
8×4 Array{Union{Missing, Char},2}:
 'a'  'a'      'a'      'a'
 'b'  'c'      'b'      'b'
 'c'  'b'      'c'      'c'
 'b'  'c'      'h'      'b'
 'c'  'f'      'i'      'c'
 'd'  'g'      'm'      'j'
 'e'  missing  'n'      'k'
 'm'  missing  missing  'm'

juila> NetMSA.matrixalignment(M)
9×4 Array{Union{Missing, Char},2}:
 'a'  'a'  'a'  'a'
 'b'  '-'  'b'  'b'
 'c'  'c'  'c'  'c'
 'b'  'b'  '-'  'b'
 'c'  'c'  '-'  'c'
 'd'  'f'  'h'  'j'
 'e'  'g'  'i'  'k'
 'm'  '-'  'm'  'm'
 '-'  '-'  'n'  '-'
```
"""
function matrixalignment(M)
  for (index, row) in enumerate(eachrow(M))
    g = rowalignment(index, M)
    if !isnothing(g)
      M = flydown(g, M, stride=g.best.row - g.pos.row)
    end
  end
  replace!(M, missing => '-')
end

end
