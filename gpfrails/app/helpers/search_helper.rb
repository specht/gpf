module SearchHelper
	def SearchHelper.text_field_with_label ak_Parameter
		ls_Result = "<table class='query-form'><tr>"
		ls_Result += "<td style='text-align: left;'><label for='#{ak_Parameter['id']}'>#{ak_Parameter['label']}:</label></td><td style='text-align: right;'><input type='text' id='#{ak_Parameter['id']}' name='#{ak_Parameter['id']}'"
		if (ak_Parameter['type'] == 'int')
			ls_Result += " size='5'"
			ls_Result += " style='text-align: right;'"
		end
		ls_Result += " value='#{ak_Parameter['default']}'"
		ls_Result += "/></td></tr></table>\n"
		return ls_Result
	end

	def text_field_with_label ak_Parameter
		return (SearchHelper.text_field_with_label ak_Parameter)
	end

	def SearchHelper.choice_field_with_label ak_Parameter
		ls_Result = "<table class='query-form'><tr>"
		ls_Result += "<td style='text-align: left;'><label for='#{ak_Parameter['id']}'>#{ak_Parameter['label']}:</label></td><td style='text-align: right;'>";
		ls_Result += "<select id='#{ak_Parameter['id']}' name='#{ak_Parameter['id']}'>"
		ak_Parameter['type'].each { |lk_Choice|
			ls_Choice = lk_Choice
			ls_Key = ls_Choice
			ls_Selected = ''
			if (lk_Choice.class == Hash)
				ls_Key = lk_Choice.keys.first
				ls_Choice = lk_Choice[ls_Key]
				ls_Choice.gsub!('Chlamydomonas', 'C.')
			end
			if (ls_Key == ak_Parameter['default'])
				ls_Select = " selected='selected'"
			end
			ls_Result += "<option value='#{ls_Key}'#{ls_Select}>#{ls_Choice}</option>"
		}
		ls_Result += "</select>"
		ls_Result += "</td></tr></table>\n"
		return ls_Result
	end
end
