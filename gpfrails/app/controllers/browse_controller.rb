class BrowseController < ApplicationController

	require 'gpf_constants'
	
	def showHit
		# take an assembly and construct a viewing window, then redirect to show
		ls_Assembly = params['assembly']
		lk_Parts = ls_Assembly.split(';')
		ls_Genome = lk_Parts[0]
		lb_Forward = lk_Parts[1][0, 1] == '+'
		lk_AssemblyParts = lk_Parts[1][1, lk_Parts[1].length - 1].split(',')
		
		lk_Result = YAML::load(Net::HTTP.get(URI.parse("#{$gs_GpfServer}?browse&genome=#{ls_Genome}&position=#{lk_AssemblyParts[0].split(':')[0]}&length=1")))
		li_ScaffoldStartPosition = lk_Result['scaffoldStart']
		li_ScaffoldLastPosition = lk_Result['scaffoldStart'] + lk_Result['scaffoldLength'] - 1

		li_Position = 0
		li_Length = $gi_BrowseWidth
		
		if (lk_AssemblyParts.size == 1)
			# immediate hit
			lk_AssemblyPart = lk_AssemblyParts[0].split(':')
			li_AssemblyPartPosition = lk_AssemblyPart[0].to_i
			li_AssemblyPartLength = lk_AssemblyPart[1].to_i
			li_Center = lb_Forward ? (li_AssemblyPartPosition + li_AssemblyPartLength / 2) : (li_AssemblyPartPosition - li_AssemblyPartLength / 2)
			
			li_Position = li_Center

			li_Position -= $gi_BrowseWidth / 2
			li_Length = $gi_BrowseWidth
		else
			# intron split hit
			lk_AssemblyPart = lk_AssemblyParts[0].split(':')
			li_AssemblyPartPosition = lk_AssemblyPart[0].to_i

			lk_AssemblyPart = lk_AssemblyParts[1].split(':')
			if (lb_Forward)
				li_AssemblyPartLength = lk_AssemblyPart[0].to_i + lk_AssemblyPart[1].to_i - li_AssemblyPartPosition
			else
				li_AssemblyPartLength = li_AssemblyPartPosition - (lk_AssemblyPart[0].to_i - lk_AssemblyPart[1].to_i)
			end
			
			li_AssemblyPartPosition -= li_AssemblyPartLength - 1 if (!lb_Forward)
			
			li_Position = li_AssemblyPartPosition - 9
			li_Length = li_AssemblyPartLength + 9
			li_Length += $gi_BrowseWidth - (li_Length % $gi_BrowseWidth)
		end
		
		li_Position -= $gi_BrowseWidth
		li_Length += $gi_BrowseWidth * 2
		
		# adjust start if necessary
		li_Position = li_ScaffoldStartPosition if (li_Position < li_ScaffoldStartPosition)
		#adjust end if necessary
		li_Length = li_ScaffoldLastPosition - li_Position + 1 if (li_Position + li_Length - 1 > li_ScaffoldLastPosition)
		
		params['genome'] = ls_Genome
		params['position'] = li_Position.to_s
		params['length'] = li_Length.to_s
		params['forward'] = lb_Forward.to_s
		
		# now redirect to show
		show()
		render :action => 'show'
	end
	
	def show
		# take a genome, a position, a length, and a direction plus an optional assembly to highlight
		
		# adjust viewing window
		li_Position = params['position'].to_i
		li_Length = params['length'].to_i
		
		lk_Result = YAML::load(Net::HTTP.get(URI.parse("#{$gs_GpfServer}?browse&genome=#{params['genome']}&position=#{li_Position}&length=1")))
		li_Position += (3 - (lk_Result['scaffoldPosition'] % 3)) % 3

		lb_Forward = params['forward'] == 'true'
		
		# fetch nucleotides from GPF server
		lk_Result = YAML::load(Net::HTTP.get(URI.parse("#{$gs_GpfServer}?browse&genome=#{params['genome']}&position=#{li_Position}&length=#{li_Length}")))
		
		# reverse and transpose DNA if necessary
		if (!lb_Forward)
			ls_Dna = lk_Result['dna'].reverse
			lk_Result['dna'] = ''
			lk_Transpose = {'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', '.' => '.'}
			(0...ls_Dna.length).each { |i| lk_Result['dna'] += lk_Transpose[ls_Dna[i, 1]] }
		end
		
		ls_Dna = lk_Result['dna']
		@mk_Result = lk_Result
		
		# translate nucleotides to amino acids
		lk_ReadingFrames = Array.new()
		(0..2).each do |li_ReadingFrame|
			ls_ReadingFrame = ''
			(0..li_ReadingFrame).each { |i| ls_ReadingFrame += '_' }
			(0...ls_Dna.length / 3).each { |i| ls_ReadingFrame += $gk_TripletToAminoAcid[ls_Dna[i * 3 + li_ReadingFrame, 3]] + '__' if (i * 3 + li_ReadingFrame + 2 < ls_Dna.length)}
			lk_ReadingFrames.push(ls_ReadingFrame)
		end
		
		@mi_Position = li_Position
		@mi_Length = li_Length
		@mk_Markup = Array.new()

		if (params['assembly'])
			ls_Assembly = params['assembly']
			lk_Parts = ls_Assembly.split(';')
			ls_Genome = lk_Parts[0]
			lb_AssemblyForward = lk_Parts[1][0, 1] == '+'
			if (lb_AssemblyForward == lb_Forward)
				@ms_PeptideShowing = params['peptide']
				@ms_AssemblyShowing = params['assembly']
				lk_AssemblyParts = lk_Parts[1][1, lk_Parts[1].length - 1].split(',')
				lk_AssemblyParts.each do |ls_AssemblyPart|
					lk_AssemblyPart = ls_AssemblyPart.split(':')
					li_PartPosition = lk_AssemblyPart[0].to_i
					li_PartLength = lk_AssemblyPart[1].to_i
					if (lb_Forward)
						li_StartOffset = li_PartPosition - @mi_Position
						li_EndOffset = li_PartPosition - @mi_Position + li_PartLength - 1
						@mk_Markup.push({'first' => li_StartOffset, 'last' => li_EndOffset, 'frame' => ls_AssemblyPart == lk_AssemblyParts.first ? li_StartOffset % 3 : (li_EndOffset + 1) % 3 })
					else
						li_StartOffset = (@mi_Position + @mi_Length - 1) - li_PartPosition
						li_EndOffset = (@mi_Position + @mi_Length - 1) - li_PartPosition + li_PartLength - 1
						@mk_Markup.push({'first' => li_StartOffset, 'last' => li_EndOffset, 'frame' => ls_AssemblyPart == lk_AssemblyParts.first ? li_StartOffset % 3 : (li_EndOffset + 1) % 3 })
					end
				end
			end
		end
		
		@mk_BrowserRows = Array.new()

		# fill @mk_BrowserRows
		(0..((ls_Dna.length - 1) / $gi_BrowseWidth)).each do |li_RowIndex|
			li_StartPosition = li_RowIndex * $gi_BrowseWidth

			lk_BrowserRow = Array.new()
			lk_BrowserRow.push(ls_Dna[li_StartPosition, $gi_BrowseWidth].dup)
			(0..2).each { |li_ReadingFrame| lk_BrowserRow.push(lk_ReadingFrames[li_ReadingFrame][li_StartPosition, $gi_BrowseWidth].dup) }

			@mk_BrowserRows.push(lk_BrowserRow)
		end
		
		# lk_Injections holds the injections for each browser row, set up
		lk_Injections = Array.new()
		(0...@mk_BrowserRows.size).each do |i|
			lk_Array = Array.new()
			(0..3).each { |i| lk_Array.push(Hash.new()) }
			lk_Injections.push(lk_Array)
		end
		
		# determine all necessary injections
		@mk_Markup.each do |lk_Markup|
			li_FirstLine = lk_Markup['first'] / $gi_BrowseWidth
			li_LastLine = lk_Markup['last'] / $gi_BrowseWidth
			(li_FirstLine .. li_LastLine).each do |li_LineIndex|
				li_StartOffset = li_LineIndex == li_FirstLine ? lk_Markup['first'] % $gi_BrowseWidth : 0
				li_EndOffset = li_LineIndex == li_LastLine ? lk_Markup['last'] % $gi_BrowseWidth + 1 : $gi_BrowseWidth
				if (li_LineIndex >= 0 && li_LineIndex < lk_Injections.size)
					# mark DNA line
					lk_Injections[li_LineIndex][0][li_StartOffset] = '<u>'
					lk_Injections[li_LineIndex][0][li_EndOffset] = '</u>'
					# mark amino acid line
					lk_Injections[li_LineIndex][lk_Markup['frame'] + 1][li_StartOffset] = '<u>'
					lk_Injections[li_LineIndex][lk_Markup['frame'] + 1][li_EndOffset] = '</u>'
				end
			end
		end
		
		# perform injections
		(0...lk_Injections.size).each do |li_RowIndex|
			(0...lk_Injections[li_RowIndex].size).each do |li_LineIndex|
				lk_Keys = lk_Injections[li_RowIndex][li_LineIndex].keys.dup
				lk_Keys.sort!
				li_TotalOffset = 0
				lk_Keys.each do |li_Offset|
					ls_Token = lk_Injections[li_RowIndex][li_LineIndex][li_Offset]
					@mk_BrowserRows[li_RowIndex][li_LineIndex].insert(li_Offset + li_TotalOffset, ls_Token)
					li_TotalOffset += ls_Token.length
				end
			end
		end
		
		# replace _ with &nbsp;
		(0...@mk_BrowserRows.size).each do |li_RowIndex|
			(0...@mk_BrowserRows[li_RowIndex].size).each do |li_LineIndex|
				@mk_BrowserRows[li_RowIndex][li_LineIndex].gsub!('_', '&nbsp;')
			end
		end
		
		@mk_Result = lk_Result
	end
end
