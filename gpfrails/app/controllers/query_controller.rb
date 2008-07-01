class QueryController < ApplicationController

	require 'gpf_constants'


	def shufflePeptide(as_Peptide)
		ls_LastChar = as_Peptide[as_Peptide.length - 1, 1]
		ls_Chopped = as_Peptide.chop
		ls_Result = ''
		while (!ls_Chopped.empty?)
			li_Index = rand(ls_Chopped.length)
			ls_Result += ls_Chopped[li_Index, 1]
			ls_Chopped[li_Index, 1] = ls_Chopped[ls_Chopped.length - 1, 1]
			ls_Chopped.chop!
		end
		ls_Result += ls_LastChar
		return ls_Result
	end


	def execute
		@ms_Result = params.to_yaml
		ls_Parameters = params['queryParameters']
		@ms_YamlResponse = Net::HTTP.get(URI.parse($gs_GpfServer + '?query=' + params['queryPeptide'] + '&' + ls_Parameters + '&fullDetails=yes'))
		lk_Result = YAML::load(@ms_YamlResponse)
		@mk_Results = lk_Result['results']
		if @mk_Results
			@mk_Results.sort! do |x, y|
				if (y['score'] == x['score'])
					if x['details']['parts'].size == y['details']['parts'].size
						if (x['peptide'] == y['peptide']) 
							if x['details']['parts'][0]['position'] == y['details']['parts'][0]['position']
								if x['details']['parts'][1]['position'] == y['details']['parts'][1]['position']
									x['details']['parts'][0]['length'] <=> y['details']['parts'][0]['length']
								else
									x['details']['parts'][1]['position'] <=> y['details']['parts'][1]['position']
								end
							else
								x['details']['parts'][0]['position'] <=> y['details']['parts'][0]['position']
							end
						else
							x['peptide'] <=> y['peptide']
						end
					else
						x['details']['parts'].size <=> y['details']['parts'].size
					end
				else
					y['score'] <=> x['score']
				end
			end
			@mk_Results.each_index do |li_Index|
				lk_Hit = @mk_Results[li_Index]
				ls_DisplayPeptide = lk_Hit['peptide'].dup
				ls_Separator = '[...]'
				if (lk_Hit['details']['parts'].size > 1)
					li_Modulo = lk_Hit['details']['parts'][0]['length'] % 3
					puts lk_Hit['assembly']
					if li_Modulo == 0
						ls_DisplayPeptide = ls_DisplayPeptide.insert(lk_Hit['details']['parts'][0]['length'] / 3, " #{ls_Separator} ")
					elsif li_Modulo == 1
						ls_DisplayPeptide = ls_DisplayPeptide.insert(lk_Hit['details']['parts'][0]['length'] / 3 + 1, ' ')
						ls_DisplayPeptide = ls_DisplayPeptide.insert(lk_Hit['details']['parts'][0]['length'] / 3, ' ' + ls_Separator)
					elsif li_Modulo == 2
						ls_DisplayPeptide = ls_DisplayPeptide.insert(lk_Hit['details']['parts'][0]['length'] / 3 + 1, ls_Separator + ' ')
						ls_DisplayPeptide = ls_DisplayPeptide.insert(lk_Hit['details']['parts'][0]['length'] / 3, ' ')
					end
				end				
				@mk_Results[li_Index]['displayPeptide'] = ls_DisplayPeptide
			end
		end

		@ms_Time = lk_Result['info']['duration']
		@ms_QueryPeptide = lk_Result['info']['query']
		@mi_MaxScore = @ms_QueryPeptide.length
		
		lk_Peptides = Hash.new()
		lk_Result['results'].each_index { |i| lk_Peptides[lk_Result['results'][i]['peptide']] = true if !lk_Peptides.has_key?(lk_Result['results'][i]['peptide']) }
		ls_Fasta = ''		
		lk_Peptides.keys.each do |ls_Peptide|
			ls_Fasta += ">#{ls_Peptide}\n#{ls_Peptide}\n"
		end
		@ms_Fasta = ls_Fasta
	end
end
