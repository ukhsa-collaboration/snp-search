
#This is the method that creates the database and imports the data.

#Called in lib/snp-search.rb

#This method guesses the reference sequence file format
def guess_sequence_format(reference_genome)
  file_extension = File.extname(reference_genome).downcase
  file_format = nil
  case file_extension
  when ".gbk", ".genbank", ".gb"
    file_format = :genbank
  when ".embl", ".emb"
    file_format = :embl
  end
  return file_format
end

# A method to populate the database with the features (genes etc) and the annotations from the gbk/embl file.  
# We include all features that are not 'source' or 'gene' as they are repetitive info.  'CDS' is the gene.
# The annotation table includes also the start and end coordinates of the CDS.  The strand is also included.  the 'locations' method is defined in bioruby under genbank.  It must be required at the top (bio).
#Also, the qualifier and value are extracted from the gbk/embl file and added to the database.
def populate_features_and_annotations(sequence_file)
  puts "Adding features and their annotations...."
  ActiveRecord::Base.transaction do
    counter = 0
    sequence_file.features.each do |feature|
      unless feature.feature == "source" || feature.feature == "gene"
        counter += 1
        puts "Total number of features and annotations added: #{counter}" if counter % 100 == 0
        db_feature = Feature.new
        db_feature.name = feature.feature
        db_feature.start = feature.locations.first.from
        db_feature.end = feature.locations.first.to
        db_feature.strand = feature.locations.first.strand
        #Add nucleotide sequence from ORIGIN of genbank file.
        db_feature.sequence = sequence_file.seq[feature.locations.first.from-1..feature.locations.first.to-1]
        db_feature.save
        # Populate the Annotation table with qualifier information from the genbank file
        feature.qualifiers.each do |qualifier|
          a = Annotation.new
          a.qualifier = qualifier.qualifier
          a.value = qualifier.value
          a.save
          db_feature.annotations << a
        end
      end
    end
  end
end

#This method populates the rest of the information, i.e. SNP information, Alleles and Genotypes.
def populate_snps_alleles_genotypes(vcf_file, cuttoff_ad)

  puts "Adding SNPs........"
  # open vcf file and parse each line
  File.open(vcf_file) do |f|
    # header names
    while line = f.gets
      if  line =~ /CHROM/
        line.chomp!
        column_headings = line.split("\t")
        strain_names = column_headings[9..-1]
        strain_names.map!{|name| name.sub(/\..*/, '')}

        strain_names.each do |str|
          ss = Strain.new
          ss.name = str
          ss.save
        end

        strains = Array.new
        strain_names.each do |strain_name|
          strain = Strain.find_by_name(strain_name) # equivalent to Strain.find.where("strains.name=?", strain_name).first
          strains << strain
        end

        good_snps = 0
        # start parsing snps
        while line = f.gets
          line.chomp!
          details = line.split("\t")
          ref = details[0]
          ref_pos = details[1].to_i

          ref_base = details[3]
          snp_bases = details[4].split(",")
          snp_qual = details [5]
          format = details[8].split(":")
          gt_array_position = format.index("GT")
          gq_array_position = format.index("GQ")
          ad_array_position = format.index("AD")
          # dp = format.index("DP")
          samples = details[9..-1]

          gts = []
          gqs = []
          ad_ratios = []

          
          next if samples.any?{|sample| sample =~ /\.\/\./}  # no coverage in at least one sample
          samples.map do |sample|
            format_values = sample.split(":") # output (e.g.): ["0/0 ", "0,255,209", "99"]
            gt = format_values[gt_array_position] # e.g.          
            gt = gt.split("/")
            next if gt.size > 1 && (gt.first != gt.last) # if its 0/1, 1/2 etc then ignore
            next if gt.first == "." # no coverage
            gt = gt.first.to_i

            gq = format_values[gq_array_position].to_f

            if ad_array_position
              # If there is AD in vcf.  Typically AD is Allele specific depth. i.e. if ref is 'A' and alt is 'G' and AD is '6,9' you got 6 A reads and 9 G reads.
              # ad below is 6 and 9 in the example above.
              ad = format_values[ad_array_position].split(",").map{|ad_value| ad_value.to_i}
              # Find the sum of all bases (sum_of_ad) reported by the ad, so its 15 in the example.
              sum_of_ad = ad.inject{|sum,x| sum + x }
              ad_ratios << ad[gt]/sum_of_ad.to_f
            end

            gqs << gq
            gts << gt
          end

          next if ad_ratios.any?{|ad_ratio| ad_ratio < cuttoff_ad.to_i} # exclude if any samples have a call ratio of less than a cuttoff set by user
          if gts.size == samples.size # if some gts have been rejected due to heterozygote or no coverage
            good_snps +=1
            
            # populate snps

            ActiveRecord::Base.transaction do
              s = Snp.new
              s.ref_pos = ref_pos
              s.qual = snp_qual
              s.save

              #  create ref allele
              ref_allele = Allele.new
              ref_allele.base = ref_base
              ref_allele.snp = s
              ref_allele.save

              s.reference_allele = ref_allele
              s.save

              snp_alleles = Array.new
              gts.uniq.select{|gt| gt > 0}.each do |gt|
                # create snp allele
                snp_allele = Allele.new
                snp_bases_index = gt - 1
                snp_allele.base = snp_bases[snp_bases_index]
                snp_allele.snp = s
                snp_allele.save
                snp_alleles << snp_allele
              end

              genos = []
              gts.each_with_index do |gt, index|
                genotype = Genotype.new
                genotype.strain = strains[index]
                #Adding the genotype quality with Genotype

                genotype.geno_qual = gqs[index]
                if gt == 0# wild type
                  genotype.allele = ref_allele
                else # snp type
                  genotype.allele = snp_alleles[gt - 1]
                end
                genos << genotype
              end

              # Using activerecord-import to speed up importing
              Genotype.import genos, :validate => false 
              puts "Total SNPs added so far: #{good_snps}" if good_snps % 100 == 0
              # puts "Total SNPs added so far: #{good_snps}"
            end
          end
        end
      end
    end
  end
  #Here we link the features to snps.
  puts "Linking features to SNPs"
  ActiveRecord::Base.transaction do
    Snp.all.each_with_index do |snp, index|
      puts "Total SNPs linked to features added so far: #{index}" if index % 100 == 0
      features = Feature.where("features.start <= ? AND features.end >= ?", snp.ref_pos, snp.ref_pos)

      unless features.empty?
        features.each do |feature|
          snp.features << feature
        end
      end
    end
  end
end
