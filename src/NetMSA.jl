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

  counts = countmap(row)
  delete!(counts, '-')
  max = findmax(counts)[1];
  c = length(row);
  if aligned(row)
    return w2 * max / c;
  else
    x = max == 1 ? 0 : max;
    return w1 * x / c;
  end
end

end
