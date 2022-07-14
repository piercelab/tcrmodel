$(document).ready(function(){
    var travfile;
    var trajfile;
    var trbvfile;
    var trbjfile;
    sele_AGeneSelection();
    sele_BGeneSelection();
    
    document.getElementById("sele_aspecies").onchange = function() {sele_AGeneSelection()};
    function sele_AGeneSelection() {
    	var x = document.getElementById("sele_aspecies").value;
    	if (x == "human"){
    	    travfile = $SCRIPT_ROOT + '/static/genemapdata/TRAV_HUMAN.json';
    	    trajfile = $SCRIPT_ROOT + '/static/genemapdata/TRAJ_HUMAN.json';
    	}
    	if (x == "mouse"){
    	    travfile = $SCRIPT_ROOT + '/static/genemapdata/TRAV_MOUSE.json';
    	    trajfile = $SCRIPT_ROOT + '/static/genemapdata/TRAJ_MOUSE.json';
    	}
    	if (x == "dog"){
    	    travfile = $SCRIPT_ROOT + '/static/genemapdata/TRAV_DOG.json';
    	    trajfile = $SCRIPT_ROOT + '/static/genemapdata/TRAJ_DOG.json';
    	}
    	if (x == "zebrafish"){
    	    travfile = $SCRIPT_ROOT + '/static/genemapdata/TRAV_ZEBRAFISH.json';
    	    trajfile = $SCRIPT_ROOT + '/static/genemapdata/TRAJ_ZEBRAFISH.json';
    	}
    	if (x == "sheep"){
    	    travfile = $SCRIPT_ROOT + '/static/genemapdata/TRAV_SHEEP.json';
    	    trajfile = $SCRIPT_ROOT + '/static/genemapdata/TRAJ_SHEEP.json';
    	}
    	
    	let travdropdown = $('#sele_trav');
    	travdropdown.empty();
    	travdropdown.append('<option selected="true" disabled value="">Choose TRAV</option>');
    	travdropdown.prop('selectedIndex', 1);
    	$.getJSON(travfile, function (data) {
    	    $.each(data, function (key, entry) {
    		travdropdown.append($('<option></option>').attr('value', entry.seq).text(entry.name));
    	    })
    	});
    	
    	let trajdropdown = $('#sele_traj');
    	trajdropdown.empty();
    	trajdropdown.append('<option selected="true" disabled value="">Choose TRAJ</option>');
    	trajdropdown.prop('selectedIndex', 1);
    	$.getJSON(trajfile, function (data) {
    	    $.each(data, function (key, entry) {
    		trajdropdown.append($('<option></option>').attr('value', entry.seq).text(entry.name));
    	    })
    	});
    }
    
    document.getElementById("sele_bspecies").onchange = function() {sele_BGeneSelection()};
    function sele_BGeneSelection() {
    	var y = document.getElementById("sele_bspecies").value;
    	if (y == "human"){
    	    trbvfile = $SCRIPT_ROOT + '/static/genemapdata/TRBV_HUMAN.json';
    	    trbjfile = $SCRIPT_ROOT + '/static/genemapdata/TRBJ_HUMAN.json';
    	}
    	if (y == "mouse"){
    	    trbvfile = $SCRIPT_ROOT + '/static/genemapdata/TRBV_MOUSE.json';
    	    trbjfile = $SCRIPT_ROOT + '/static/genemapdata/TRBJ_MOUSE.json';
    	}
    	if (y == "dog"){
    	    trbvfile = $SCRIPT_ROOT + '/static/genemapdata/TRBV_DOG.json';
    	    trbjfile = $SCRIPT_ROOT + '/static/genemapdata/TRBJ_DOG.json';
    	}
    	if (y == "trout"){
    	    trbvfile = $SCRIPT_ROOT + '/static/genemapdata/TRBV_RAINBOW_TROUT.json';
    	    trbjfile = $SCRIPT_ROOT + '/static/genemapdata/TRBJ_RAINBOW_TROUT.json';
    	}
    	if (y == "monkey"){
    	    trbvfile = $SCRIPT_ROOT + '/static/genemapdata/TRBV_RHESUS_MONKEY.json';
    	    trbjfile = $SCRIPT_ROOT + '/static/genemapdata/TRBJ_RHESUS_MONKEY.json';
    	}
    	
    	let trbvdropdown = $('#sele_trbv');
    	trbvdropdown.empty();
    	trbvdropdown.append('<option selected="true" disabled value="">Choose TRBV</option>');
    	trbvdropdown.prop('selectedIndex', 1);
    	$.getJSON(trbvfile, function (data) {
    	    $.each(data, function (key, entry) {
    		trbvdropdown.append($('<option></option>').attr('value', entry.seq).text(entry.name));
    	    })
    	});
    	
    	let trbjdropdown = $('#sele_trbj');
    	trbjdropdown.empty();
    	trbjdropdown.append('<option selected="true" disabled value="">Choose TRBJ</option>');
    	trbjdropdown.prop('selectedIndex', 1);
    	$.getJSON(trbjfile, function (data) {
    	    $.each(data, function (key, entry) {
    		trbjdropdown.append($('<option></option>').attr('value', entry.seq).text(entry.name));
    	    })
    	});
    }
    
    
    $("#sele_loadbtntcra").click(function () {
        var a_trav = document.getElementById("sele_trav").value;
        var a_traj = document.getElementById("sele_traj").value;
        var a_cdr3 = (document.getElementById("sele_tacdrseq").value).toUpperCase();
        var aseq = a_trav+a_cdr3+a_traj;
        document.getElementById("sele_taseq").value = aseq;
        document.getElementById("sele_trav_id").innerHTML = "";
        if( !a_cdr3 ) {
    	    document.getElementById("sele_trav_id").innerHTML = "&#9888; CDR3 sequence not entered!";
        }
        if( !a_traj ) {
    	    document.getElementById("sele_trav_id").innerHTML = "&#9888; TRAJ not selected!";
        }
        if( !a_trav ) {
    	    document.getElementById("sele_trav_id").innerHTML = "&#9888; TRAV not selected!";
        }
    });
    
    $("#sele_loadbtntcrb").click(function () {
        var b_trbv = document.getElementById("sele_trbv").value;
        var b_trbj = document.getElementById("sele_trbj").value;
        var b_cdr3 = (document.getElementById("sele_tbcdrseq").value).toUpperCase();
        var bseq = b_trbv+b_cdr3+b_trbj;
        document.getElementById("sele_tbseq").value = bseq;
        document.getElementById("sele_trbv_id").innerHTML = "";
        if( !b_cdr3 ) {
    	    document.getElementById("sele_trbv_id").innerHTML = "&#9888; CDR3 sequence not entered!";
        }
        if( !b_trbj ) {
    	    document.getElementById("sele_trbv_id").innerHTML = "&#9888; TRBJ not selected!";
        }
        if( !b_trbv ) {
    	    document.getElementById("sele_trbv_id").innerHTML = "&#9888; TRBV not selected!";
        }
    });
    
    
    $("#btn1").click(function () {
    	$("#achain").val("KTTQPISMDSYEGQEVNITCSHNNIATNDYITWYQQFPSQGPRFIIQGYKTKVTNEVASLFIPADRKSSTLSLPRVSLSDTAVYYCLVGDMDQAGTALIFGKGTTLSVS");
    	$("#bchain").val("VTQSPTHLIKTRGQQVTLRCSPKSGHDTVSWYQQALGQGPQFIFQYYEEEERQRGNFPDRFSGHQFPNYSSELNVNALLLGDSALYLCASSLGQTNYGYTFGSGTRLTVV");
    	$("#achain").keyup();
    	$("#bchain").keyup();
    });
    
    $("#btn2").click(function () {
    	$("#aspecies").prop('selectedIndex', 1);
    	$("#bspecies").prop('selectedIndex', 1);
    	$("#sele_tacdrseq").val("CAVGGSQGNLIF");
    	$("#sele_tbcdrseq").val("CASSIRSSYEQYF");
        $('select#sele_trav option:contains("TRAV8-6*01")').prop('selected',true);
        $('select#sele_traj option:contains("TRAJ42*01")').prop('selected',true);
        $('select#sele_trbv option:contains("TRBV19*01")').prop('selected',true);
        $('select#sele_trbj option:contains("TRBJ2-7*01")').prop('selected',true);
    	$("#sele_loadbtntcra").click();
        $("#sele_loadbtntcrb").click();
    });
    
    
    $("#achain").keyup(function() { 
        var x = document.getElementById("achain").value;
        x = x.replace(/(\r\n|\n|\s|\r)/gm,"");
        var found = 0;
        travfile = $SCRIPT_ROOT + '/static/genemapdata/TRAV_HUMAN.json';
        $.getJSON(travfile, function (data) {
            $.each(data, function (key, entry) {
                var re = new RegExp(entry.seq.substring(3))
        	if (re.test(x)) {
            	    document.getElementById("trav_id").innerHTML = "Human " + entry.name + " gene identified";	
            	    found = 1;
            	    return true;
        	}
            })
        });		
        if (found == 0) { document.getElementById("trav_id").innerHTML = ""; }
    });
    
    $("#bchain").keyup(function() {
        var x = document.getElementById("bchain").value;
        x = x.replace(/(\r\n|\n|\s|\r)/gm,"");
    	var found = 0;
        travfile = $SCRIPT_ROOT + '/static/genemapdata/TRBV_HUMAN.json';
        $.getJSON(trbvfile, function (data) {
            $.each(data, function (key, entry) {
                var re = new RegExp(entry.seq.substring(3))
                if (re.test(x)) {
                    document.getElementById("trbv_id").innerHTML = "Human " + entry.name + " gene identified";
                    found = 1;
                    return true;
                }
            })
        });
        if (found == 0) { document.getElementById("trbv_id").innerHTML = ""; }
    });


    $("#reset3").click(function () {	
    	document.getElementById("trav_id").innerHTML = "";
    	document.getElementById("trbv_id").innerHTML = "";
    });
});
