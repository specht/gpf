# Copyright (c) 2007-2008 Michael Specht
# 
# This file is part of GPF.
# 
# GPF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# GPF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GPF.  If not, see <http://www.gnu.org/licenses/>.

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
