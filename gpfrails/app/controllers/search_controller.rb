class SearchController < ApplicationController

	require 'search_helper'
	require 'gpf_constants'
	
	def initialize
		@mk_DefaultQueries = ['WLQYSEVIHAR', 'ALEGELYDTFK',
			'LYDEELQAIAK', 'LEEGEMPLNTYSNK', 'LAPYSEVFGLAR']
		@mb_ShowSidebar = true
		@ms_PageTitle = 'Search'
	end

	def get_Parameters
		ls_Response = Net::HTTP.get(URI.parse($gs_GpfServer + '/?getParameters'))
		lk_Parameters = YAML::load(ls_Response)
		@parameters = ''
		lk_Parameters['parameters'].each do |lk_Parameter|
			if (params[lk_Parameter['id']])
				lk_Parameter['default'] = params[lk_Parameter['id']]
			end
			if (lk_Parameter['type'] == 'int')
				@parameters += SearchHelper.text_field_with_label lk_Parameter
			elsif (lk_Parameter['type'].class == Array)
				@parameters += SearchHelper.choice_field_with_label lk_Parameter
			end
		end
	rescue Exception:
		ShowError "Unable to connect to GPF server", "<p>We are sorry for the inconvenience, but the GPF server is temporarily down.</p><p>Please try again later.</p>"
	end
	
	def index
		get_Parameters
	end
	
	def batch
		@mb_BatchProcessing = true
		@ms_PageTitle = 'Batch processing'
		get_Parameters
	end

	def query
		get_Parameters
		
		# sanitize query peptide
		if (params['peptide'] != nil)
			params['peptide'].gsub!(' ', '')
			params['peptide'].gsub!("\t", '')
			params['peptide'].gsub!("\n", '')
			params['peptide'].gsub!("\r", '')
			params['peptide'].upcase!
			@ms_Peptide = params['peptide']
		end

		@ms_QueryParameters = ''
		request.parameters.each_pair do |ls_Key, ls_Value|
			@ms_QueryParameters += '&' if !@ms_QueryParameters.empty?
			@ms_QueryParameters += "#{ls_Key}=#{ls_Value}" if ls_Value != nil && !ls_Value.empty?
		end
	end
	
	def ShowError as_Heading, as_ErrorMessage
		@ms_Heading = as_Heading
		@ms_ErrorMessage = as_ErrorMessage
		@mb_ShowSidebar = false
		@ms_PageTitle = 'Error'
		render :action => 'error'
	end
	
	def resultDownload
		headers['Content-Type'] = 'text/plain'
		headers['Content-Disposition'] = "attachment; filename=#{params['filename']}"
		@ms_Contents = params['contents']
		render :action => 'resultDownload', :layout => false
	end
	
	def error
		@ms_PageTitle = 'Error'
	end
end
