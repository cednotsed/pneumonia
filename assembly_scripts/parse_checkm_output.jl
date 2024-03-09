using JSON

# Specify the path to your file
file_path = "/mnt/c/git_repos/pneumonia/results/assembly_out/checkm_out/storage/bin_stats_ext.tsv"

# Read lines from the file
try
    lines = readlines(file_path)
    
    header = ""

    # Display the lines
    println("Lines read from the file:")
    for line in lines[1:2]
        # Splitting the data using space as the delimiter
        line = replace(line, r"\s+" => " ")
        line = replace(line, "'" => "")
        line = replace(line, r"[{}]" => "")
        split_data = split(line, " ")
        id = split_data[1]
        second_col = join(split_data[2:end])
        second_col = strip(second_col, ' ')
		
		# Parse second column
        pairs = split(second_col, ',')[1:27]
        values = [split(strip(kv), ':')[2] for kv in pairs]
        
		# Combine keys and values into TSV format
        # Convert list to TSV format
        tsv_data = join(string.(values), '\t')

        println("id: $id")
        println("$tsv_data")
    end
catch err
    println("Error reading the file: $err")
end

## Splitting the data using space as the delimiter
#split_data = split(data, ' ')
#
## Extracting information using split data
#barcode = split_data[1]
#info_dict = Dict(split(kv, ':') for kv in split_data[2:end])
#
## Displaying the results
#println("Barcode: $barcode")
#println("Info:")
#println(info_dict)
