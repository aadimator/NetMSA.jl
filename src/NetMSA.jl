module NetMSA

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

end
