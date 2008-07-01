require 'yaml'
lk_Document = YAML::load(File.open('res/GpfParameters.yaml'))

lk_File = File.open('generated/GpfParameters.h', 'w');
lk_File << "#pragma once\n\n\n"

# dump GPF parameter names
lk_File << "// GPF parameter names\n"
lk_File << "struct r_GpfParameterName\n{\n"
lk_File << "\tenum Enumeration\n\t{\n"
lk_Document['parameters'].each {
	|lk_Parameter|
	ls_Id = lk_Parameter['id']
	ls_Id[0, 1] = ls_Id[0, 1].upcase
	ls_Line = "\t\t" + ls_Id
	if (lk_Parameter != lk_Document['parameters'].last)
		ls_Line += ','
	end
	ls_Line += "\n"
	lk_File << ls_Line
}
lk_File << "\t};\n"
lk_File << "};\n"

lk_File << "\n\n";

#dump enum definitions
lk_File << "// GPF enum definitions\n"
lk_Document['enums'].each {
	|lk_Enum|
	ls_Id = lk_Enum[0]
	ls_Choices_ = lk_Enum[1]
	lk_File << "struct r_" << ls_Id << "\n"
	lk_File << "{\n"
	lk_File << "\tenum Enumeration\n\t{\n"
	ls_Choices_.each {
		|lk_Choice|
		ls_Choice = "yay"
		if (lk_Choice.class != String)
			ls_Choice = lk_Choice.keys[0].dup
		else
			ls_Choice = lk_Choice
		end
		ls_Choice[0, 1] = ls_Choice[0, 1].upcase
		lk_File << "\t\t" << ls_Choice
		if (lk_Choice != ls_Choices_.last)
			lk_File << ","
		end
		lk_File << "\n"
	}
	lk_File << "\t};\n};\n\n\n"
}

lk_File << "// GPF enum options\n"
lk_Document['enums'].each {
# QMap<QString, r_YesNoOption::Enumeration> gk_YesNoOptions;
	|lk_Enum|
	ls_Id = lk_Enum[0]
	ls_Choices_ = lk_Enum[1]
	lk_File << "extern QHash<QString, int> gk_" << ls_Id << "s;\n"
	lk_File << "extern QHash<int, QString> gk_" << ls_Id << "sReverse;\n"
}

lk_File.close()

lk_File = File.open('generated/GpfOptionsDeclare.inc.cpp', 'w')

lk_File << "// GPF enum options\n"
lk_Document['enums'].each {
# QMap<QString, r_YesNoOption::Enumeration> gk_YesNoOptions;
	|lk_Enum|
	ls_Id = lk_Enum[0]
	ls_Choices_ = lk_Enum[1]
	lk_File << "QHash<QString, int> gk_" << ls_Id << "s;\n"
	lk_File << "QHash<int, QString> gk_" << ls_Id << "sReverse;\n"
}

lk_File.close()


lk_File = File.open('generated/GpfParametersInitialize.inc.cpp', 'w')

# mk_Parameters[r_GpfParameterName::PrecursorMass] = new k_GpfParameterInt(id, label, description, default);
lk_File << "// GPF parameter intializations\n"
lk_Document['parameters'].each {
	|lk_Parameter|
	ls_Id = lk_Parameter['id']
	ls_Id[0, 1] = ls_Id[0, 1].upcase
	ls_Variable = "mk_" + ls_Id + "Parameter"
	if (lk_Parameter['description'] == nil)
		lk_Parameter['description'] = ""
	end
	
	ls_Line = "mk_Parameters[r_GpfParameterName::" + ls_Id + "] = new k_GpfParameter";
	lk_Parameter['id'][0, 1] = lk_Parameter['id'][0, 1].downcase
	
	lb_Enum = false

	case lk_Parameter['type']
	when "float"
		ls_Line += "Double"
	when "string"
		ls_Line += "String"
	when "int"
		ls_Line += "Int"
	else 
		ls_Line += "Enum"
		lb_Enum = true
	end
	
	ls_Line += "(\"" + ls_Id + "\", \"" + lk_Parameter['label'] + "\", \"" + lk_Parameter['description'] + "\", "

	if (lk_Parameter['type'] == "string")
		ls_Line += "\""
	end
	
	ls_EnumName = ""
	if (lk_Parameter['type'].class == Array)
		ls_Line += "r_" 
		lk_Document['enums'].each_pair {
			|lk_Key, lk_Value|
			if (lk_Value == lk_Parameter['type'])
				ls_Line += lk_Key
				ls_EnumName = lk_Key
			end
		}
		ls_Line += "::"
	end
	
	ls_Value = lk_Parameter['default'].to_s.dup
	ls_Value[0, 1] = ls_Value[0, 1].upcase
	ls_Line += ls_Value
	
	if (lk_Parameter['type'] == "string")
		ls_Line += "\""
	end
	
	if (lb_Enum)
		ls_Line += ", gk_" + ls_EnumName + "s, gk_" + ls_EnumName + "sReverse"
	end

	ls_Line += ");\n"
	lk_File << ls_Line
}

lk_File << "\n// fill GPF enum option arrays\n"
lk_Document['enums'].each {
	|lk_Enum|
	ls_Id = lk_Enum[0]
	ls_Choices_ = lk_Enum[1]
	li_Counter = 0;
	ls_Choices_.each {
		|lk_Choice|
		ls_Key = lk_Choice.dup
		if (ls_Key.class == Hash)
			ls_Key = ls_Key.keys[0].dup
		end
		ls_Key[0, 1] = ls_Key[0, 1].downcase
		lk_File << "gk_" << ls_Id << "s[\"" << ls_Key << "\"] = " << li_Counter << ";\n"
		lk_File << "gk_" << ls_Id << "sReverse[" << li_Counter << "] = \"" << ls_Key << "\";\n"
		li_Counter = li_Counter + 1
	}
}

lk_File.close();
