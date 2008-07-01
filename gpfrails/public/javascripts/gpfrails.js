function toggleMore(ai_Id)
{
	var lb_Explorer = BrowserDetect.browser == "Explorer";
	if (lb_Explorer)
	{
		lk_Elements = document.getElementsByName(ai_Id);
		var ls_NewStyle = lk_Elements[0].style.getAttribute("display") == "none"? "inline": "none";
		for (var i = 0; i < lk_Elements.length; ++i)
			lk_Elements[i].style.setAttribute("display", ls_NewStyle, true);
	}
	else
	{
		lk_Elements = document.getElementsByName(ai_Id);
		var ls_NewStyle = lk_Elements[0].style.display == "none"? "table-row": "none";
		for (var i = 0; i < lk_Elements.length; ++i)
			lk_Elements[i].style.display = ls_NewStyle;
	}
}
