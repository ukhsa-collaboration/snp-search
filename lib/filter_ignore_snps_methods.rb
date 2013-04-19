
# require "/Volumes/NGS2_DataRAID/projects/ali/GAS/snp-search-2.0.0/lib/information_methods.rb"
# This method performs several queries to ignore elements of the data for fasta or tabular output.
# Its is called in lib/snp-search.rb

require 'output_information_methods'

def get_snps(out, ignore_snps_on_annotation, ignore_snps_in_range, ignore_strains, remove_non_informative_snps, fasta_output, tabular_output, cuttoff_genotype, cuttoff_snp, tree, fasttree_path)

  strains = Strain.all

  sequence_hash = Hash.new
  sequence_hash["ref"] = Array.new
  strains.each do |strain|
    sequence_hash[strain.name] = Array.new
  end

  snps_array = Array.new
  snp_positions = Array.new

  # output opened for data input
  output = File.open(out, "w")
  tab_delim_file_name = File.basename(out, File.extname(out)) + "_snps.tsv"
  tab_delim_file = File.open(tab_delim_file_name, "w")
  position_map_file_name = File.basename(out, File.extname(out)) + "_snps_positions.txt"
  position_map_file = File.open(position_map_file_name, "w")


  snps_within_features_with_annotation = ""
  # Perform query
  # puts ignore_snps_on_annotation.inspect
  if ignore_snps_on_annotation
    annotations_where = ignore_snps_on_annotation.split(",").map{|annotation| "annotations.value LIKE '%#{annotation}%'"}.join(" OR ")
    features_with_annotation = Feature.joins(:annotations).where(annotations_where)
    snps_within_features_with_annotation = Snp.joins(:features).where("features.id IN (?)", features_with_annotation.collect{|feature| feature.id})
  end

  if snps_within_features_with_annotation.empty?
    snps = Snp.all
  else
    snps = Snp.where("snps.id NOT IN (?)", snps_within_features_with_annotation.collect{|snp| snp.id})
  end

  positions_to_ignore = Array.new
  if ignore_snps_in_range
    range_strings = ignore_snps_in_range.split(",")
    range_strings.each do |range|
      start_position, end_position = range.split("..")
      positions_to_ignore += (start_position.to_i..end_position.to_i).to_a
    end
  end

  if ignore_strains
    strains_to_ignore = ignore_strains.split(",")
  end


  i = 0
  puts "Your Query is submitted and is being processed......."
  strains = Strain.find(:all)
  if ignore_strains
    strains_to_ignore = ignore_strains.split(",")
    strains.reject!{|strain| strains_to_ignore.include?(strain.name)}
  end

  snps.each do |snp|

    ActiveRecord::Base.transaction do
      i += 1
      next if positions_to_ignore.include?(snp.ref_pos) # Ignore positions that user specified
      alleles = snp.alleles
            
      genotypes = snp.alleles.collect{|allele| allele.genotypes}.flatten

      snp_qual = Snp.find_by_sql("select qual from snps where snps.id = #{snp.id}")
      # ignore snp if the snp qual is less than cuttoff.
      next if snp_qual.any?{|snps_quality| snps_quality.qual < cuttoff_snp.to_i}
      
      next if alleles.any?{|allele| allele.base.length > 1} # indel
      next unless genotypes.all?{|genotype| genotype.geno_qual >= cuttoff_genotype} # all geno quals > cutoff
      # puts "#{i} SNPs processed so far" if i % 100 == 0
      strain_alleles = Hash.new
      strains.each do |strain|
        strain_genotype = genotypes.select{|genotype| genotype.strain_id == strain.id}.first
        strain_allele = alleles.select{|allele| allele.id == strain_genotype.allele_id}.first

        strain_alleles[strain.name] = strain_allele.base
     end

      if remove_non_informative_snps
        next if strain_alleles.values.uniq.size == 1 # remove non-informative SNPs
      end

      snp_positions << snp.ref_pos
      snps_array << snp
      strain_alleles.each do |strain_name, allele_base|
        sequence_hash[strain_name] << allele_base
      end
      sequence_hash["ref"] << snp.reference_allele.base
    end
  end

  # If user has specified a tabular output
  if tabular_output
    output_information_methods(snps_array, output, cuttoff_genotype, cuttoff_snp, true)
  # If user has specified a fasta output
  elsif fasta_output
    # generate FASTA file
    output.puts ">ref\n#{sequence_hash["ref"].join("")}"
    tab_delim_file.puts "\t#{snp_positions.join("\t")}"
    tab_delim_file.puts "ref\t#{sequence_hash["ref"].join("\t")}"
    strains.each do |strain|
      output.puts ">#{strain.name}\n#{sequence_hash[strain.name].join("")}"
      tab_delim_file.puts "#{strain.name}\t#{sequence_hash[strain.name].join("\t")}"
    end

    snp_positions.each_with_index do |snp_position, index|
      position_map_file.puts "#{index+1} => #{snp_position}"
    end
  end
  # If user has chosen a newick output.
  if tree
    nwk_out_file_name = File.basename(out, File.extname(out)) + ".nwk"
    puts "running phylogeny"
    `#{fasttree_path} -fastest -nt #{output} > #{nwk_out_file_name}`
  end

  output.close
  tab_delim_file.close
end
