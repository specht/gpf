class HelpController < ApplicationController

	def index
		@ms_PageTitle = 'Help'
	end
	
	def faq
		@ms_PageTitle = 'Frequently asked questions'
	end
	
	def introduction
		@ms_PageTitle = 'Introduction to peptide searching'
	end
	
end
