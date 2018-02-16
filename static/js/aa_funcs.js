$(document).ready(function() {
	$('.ali').each(function(index) {
                //var $this = $(this);
                //$this.empty();
		var s1 = document.getElementsByClassName("t1ali")[index].innerHTML;
		var s2 = document.getElementsByClassName("t2ali")[index].innerHTML;
		var s3 = "";
		
		for(var i = 0; i < s1.length; i++) {
		    var c1 = s1.charAt(i);
		    var c2 = s2.charAt(i);
		    //REF: https://www.ebi.ac.uk/Tools/msa/clustalo/help/faq.html
		    //asterisk
		    if (c1 == c2) { s3 = s3 + "*"; }
		    //colon
		    else if ( ("STA".indexOf(c1) >= 0) &&("STA".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
		    else if ( ("NEQK".indexOf(c1) >= 0) &&("NEQK".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
                    else if ( ("NHQK".indexOf(c1) >= 0) &&("NHQK".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
                    else if ( ("NDEQ".indexOf(c1) >= 0) &&("NDEQ".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
                    else if ( ("QHRK".indexOf(c1) >= 0) &&("QHRK".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
                    else if ( ("MILV".indexOf(c1) >= 0) &&("MILV".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
                    else if ( ("MILF".indexOf(c1) >= 0) &&("MILF".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
                    else if ( ("HY".indexOf(c1) >= 0) &&("HY".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
                    else if ( ("FYW".indexOf(c1) >= 0) &&("FYW".indexOf(c2) >= 0) ) { s3 = s3 + ":"; }
		    //period
                    else if ( ("CSA".indexOf(c1) >= 0) &&("CSA".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("ATV".indexOf(c1) >= 0) &&("ATV".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("SAG".indexOf(c1) >= 0) &&("SAG".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("STNK".indexOf(c1) >= 0) &&("STNK".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("STPA".indexOf(c1) >= 0) &&("STPA".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("SGND".indexOf(c1) >= 0) &&("SGND".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("SNDEQK".indexOf(c1) >= 0) &&("SNDEQK".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("NDEQHK".indexOf(c1) >= 0) &&("NDEQHK".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("NEQHRK".indexOf(c1) >= 0) &&("NEQHRK".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("FVLIM".indexOf(c1) >= 0) &&("FVLIM".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
                    else if ( ("HFY".indexOf(c1) >= 0) &&("HFY".indexOf(c2) >= 0) ) { s3 = s3 + "."; }
		    //space
		    else { s3 = s3 + "&nbsp;"; }
		}
		$(this).find(".ast").html(s3);
	    });
	
	$('.aaseqali').each(function (index) {
		var characters = $(this).text().split("");
		$this = $(this);
		$this.empty();
		$.each(characters, function (i, el) {
			$this.append("<span class=\""+el+"\">" + el + "</span");
		    });
	    });
	
	
    });

