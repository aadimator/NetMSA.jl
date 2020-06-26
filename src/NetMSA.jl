module NetMSA

using StatsBase

# Write your package code here.
function createPeerMatrix(inputStrings::Array{String,1} ):: Matrix{Union{Missing, Char}}
		col_size = size(inputStrings, 1);
		row_size = max([length(s) for s in inputStrings]...);
		peer_matrix = Matrix{Union{Missing, Char}}(missing, row_size, col_size);
		for (i, s) in enumerate(inputStrings)
			peer_matrix[1: length(s), i] = collect(s);
		end
		return peer_matrix;
end


function mostfrequent(row)
  counts = countmap(row)
  delete!(counts, '-')
  max = findmax(counts);
  return max;
end

function aligned(row)::Bool
  row = Set(row)
  return length(row) == 1 || (length(row) == 2 && ('-' in row || missing in row))
end

function full(row)::Bool
  return length(Set(row)) == 1
end

function weight(row, w1=0.25, w2=0.5, w3=1.0)
  if full(row)
    return w3;
  end

  max = mostfrequent(row)[1];
  c = length(row);
  if aligned(row)
    return w2 * max / c;
  else
    x = max == 1 ? 0 : max;
    return w1 * x / c;
  end
end

function objective(M, rowind::Int64; endind::Int64=-1)
  weights = sum(weight.(eachrow(M[rowind:end, :])))
  C = mostfrequent(M[rowind, :])[1];
  A = sum(aligned.(eachrow(M))[rowind:end])

  endind = endind == -1 ? size(M)[1] : endind;
  if endind > size(M)[1]
    throw(ArgumentError("endind exceeds the matrix size"));
  end
  counts = countmap(M[rowind:endind, :]);
  Gaps = get(counts, '-', 0);

  return weights * (A * C)/(1 + Gaps)
end

end
