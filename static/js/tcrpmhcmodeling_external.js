$(document).ready(function(){
    var travfile;
    var trajfile;
    var trbvfile;
    var trbjfile;
    AGeneSelection();
    BGeneSelection();
    sele_AGeneSelection();
    sele_BGeneSelection();

    Resetforms();
    function Resetforms() {
	$("#reset1").click();
    }

    document.getElementById("aspecies").onchange = function() {AGeneSelection()};
    function AGeneSelection() {
	var x = document.getElementById("aspecies").value;
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
	
	let travdropdown = $('#trav');
	travdropdown.empty();
	travdropdown.append('<option selected="true" disabled value="None">Choose TRAV</option>');
	travdropdown.prop('selectedIndex', 1);
	$.getJSON(travfile, function (data) {
	    $.each(data, function (key, entry) {
		travdropdown.append($('<option></option>').attr('value', entry.seq).text(entry.name));
	    })
		});
	
	let trajdropdown = $('#traj');
	trajdropdown.empty();
	trajdropdown.append('<option selected="true" disabled value="">Choose TRAJ</option>');
	trajdropdown.prop('selectedIndex', 1);
	$.getJSON(trajfile, function (data) {
	    $.each(data, function (key, entry) {
		trajdropdown.append($('<option></option>').attr('value', entry.seq).text(entry.name));
	    })
		});
    }

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
	travdropdown.append('<option selected="true" disabled value="None">Choose TRAV</option>');
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
    
    document.getElementById("bspecies").onchange = function() {BGeneSelection()};
    document.getElementById("sele_bspecies").onchange = function() {sele_BGeneSelection()};
    function BGeneSelection() {
	var y = document.getElementById("bspecies").value;
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
	
	let trbvdropdown = $('#trbv');
	trbvdropdown.empty();
	trbvdropdown.append('<option selected="true" disabled value="None">Choose TRBV</option>');
	trbvdropdown.prop('selectedIndex', 1);
	$.getJSON(trbvfile, function (data) {
	    $.each(data, function (key, entry) {
		trbvdropdown.append($('<option></option>').attr('value', entry.seq).text(entry.name));
	    })
		});
	
	let trbjdropdown = $('#trbj');
	trbjdropdown.empty();
	trbjdropdown.append('<option selected="true" disabled value="">Choose TRBJ</option>');
	trbjdropdown.prop('selectedIndex', 1);
	$.getJSON(trbjfile, function (data) {
	    $.each(data, function (key, entry) {
		trbjdropdown.append($('<option></option>').attr('value', entry.seq).text(entry.name));
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
	trbvdropdown.append('<option selected="true" disabled value="None">Choose TRBV</option>');
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
    

    var mhc1afile;
    document.getElementById("mhc1speciestype").onchange = function() {MHC1SpeciesSelection()};
    MHC1SpeciesSelection();
    function MHC1SpeciesSelection() {
	var x = document.getElementById("mhc1speciestype").value;
	if (x == "human"){
            mhc1afile = $SCRIPT_ROOT + '/static/genemapdata/hla_ref_set_class_i.json';
	}
	let mhc1adropdown = $('#mhc1a');
	mhc1adropdown.empty();
	mhc1adropdown.append('<option selected="true" disabled value="">Choose MHC class I</option>');
	mhc1adropdown.prop('selectedIndex', 0);
	$.getJSON(mhc1afile, function (data) {
	    $.each(data, function (key, entry) {
		mhc1adropdown.append($('<option></option>').attr('value', entry.fullseq).text(entry.ref_name1));
	    })
		});
    }

    document.getElementById("mhc1a").onchange = function() {MHC1GeneSelection()};
    MHC1GeneSelection();
    function MHC1GeneSelection() {
	var mhc1aseq = document.getElementById("mhc1a").value;
	document.getElementById("m1aseq").value = mhc1aseq;
	$("#m1aseq").val(mhc1aseq);
    }

    document.getElementById("F1mhc1speciestype").onchange = function() {F1MHC1SpeciesSelection()};
    F1MHC1SpeciesSelection();
    function F1MHC1SpeciesSelection() {
	var mhc1afile;
	var x = document.getElementById("F1mhc1speciestype").value;
	if (x == "human"){
            mhc1afile = $SCRIPT_ROOT + '/static/genemapdata/hla_ref_set_class_i.json';
	}
	let mhc1adropdown = $('#F1mhc1a');
	mhc1adropdown.empty();
	mhc1adropdown.append('<option selected="true" disabled value="">Choose MHC class I</option>');
	mhc1adropdown.prop('selectedIndex', 0);
	$.getJSON(mhc1afile, function (data) {
	    $.each(data, function (key, entry) {
		mhc1adropdown.append($('<option></option>').attr('value', entry.fullseq).text(entry.ref_name1));
	    })
		});
    }

    document.getElementById("F1mhc1a").onchange = function() {F1MHC1GeneSelection()};
    F1MHC1GeneSelection();
    function F1MHC1GeneSelection() {
	var mhc1aseq = document.getElementById("F1mhc1a").value;
	document.getElementById("F1m1aseq").value = mhc1aseq;
	$("#F1m1aseq").val(mhc1aseq);
    }


    var mhc2afile;
    var mhc2bfile;
    document.getElementById("mhc2speciestype").onchange = function() {MHC2SpeciesSelection()};
    MHC2SpeciesSelection();
    function MHC2SpeciesSelection() {
	var x = document.getElementById("mhc2speciestype").value;
	if (x == "human"){
            mhc2afile = $SCRIPT_ROOT + '/static/genemapdata/hla_ref_set_class_iiA.json';
            mhc2bfile = $SCRIPT_ROOT + '/static/genemapdata/hla_ref_set_class_iiB.json';
	}
	let mhc2adropdown = $('#mhc2a');
	mhc2adropdown.empty();
	mhc2adropdown.append('<option selected="true" disabled value="">Choose MHC II ??</option>');
	mhc2adropdown.prop('selectedIndex', 0);
	$.getJSON(mhc2afile, function (data) {
	    $.each(data, function (key, entry) {
		mhc2adropdown.append($('<option></option>').attr('value', entry.fullseq).text(entry.ref_name1));
	    })
		});
	let mhc2bdropdown = $('#mhc2b');
	mhc2bdropdown.empty();
	mhc2bdropdown.append('<option selected="true" disabled value="">Choose MHC II ??</option>');
	mhc2bdropdown.prop('selectedIndex', 0);
	$.getJSON(mhc2bfile, function (data) {
	    $.each(data, function (key, entry) {
		mhc2bdropdown.append($('<option></option>').attr('value', entry.fullseq).text(entry.ref_name1));
	    })
		});
    }    
    document.getElementById("mhc2a").onchange = function() {MHC2GeneSelection()};
    document.getElementById("mhc2b").onchange = function() {MHC2GeneSelection()};
    MHC2GeneSelection();
    function MHC2GeneSelection() {
	var mhc2aseq = document.getElementById("mhc2a").value;
	document.getElementById("m2aseq").value = mhc2aseq;
	$("#m2aseq").val(mhc2aseq);
	var mhc2bseq = document.getElementById("mhc2b").value;
	document.getElementById("m2bseq").value = mhc2bseq;
	$("#m2bseq").val(mhc2bseq);
    }

    document.getElementById("F1mhc2speciestype").onchange = function() {F1MHC2SpeciesSelection()};
    F1MHC2SpeciesSelection();
    function F1MHC2SpeciesSelection() {
	var mhc2afile;
	var mhc2bfile;
	var x = document.getElementById("F1mhc2speciestype").value;
	if (x == "human"){
            mhc2afile = $SCRIPT_ROOT + '/static/genemapdata/hla_ref_set_class_iiA.json';
            mhc2bfile = $SCRIPT_ROOT + '/static/genemapdata/hla_ref_set_class_iiB.json';
	}
	let mhc2adropdown = $('#F1mhc2a');
	mhc2adropdown.empty();
	mhc2adropdown.append('<option selected="true" disabled value="">Choose MHC II ??</option>');
	mhc2adropdown.prop('selectedIndex', 0);
	$.getJSON(mhc2afile, function (data) {
	    $.each(data, function (key, entry) {
		mhc2adropdown.append($('<option></option>').attr('value', entry.fullseq).text(entry.ref_name1));
	    })
		});
	let mhc2bdropdown = $('#F1mhc2b');
	mhc2bdropdown.empty();
	mhc2bdropdown.append('<option selected="true" disabled value="">Choose MHC II ??</option>');
	mhc2bdropdown.prop('selectedIndex', 0);
	$.getJSON(mhc2bfile, function (data) {
	    $.each(data, function (key, entry) {
		mhc2bdropdown.append($('<option></option>').attr('value', entry.fullseq).text(entry.ref_name1));
	    })
		});
    }    
    document.getElementById("F1mhc2a").onchange = function() {F1MHC2GeneSelection()};
    document.getElementById("F1mhc2b").onchange = function() {F1MHC2GeneSelection()};
    F1MHC2GeneSelection();
    function F1MHC2GeneSelection() {
	var mhc2aseq = document.getElementById("F1mhc2a").value;
	document.getElementById("F1m2aseq").value = mhc2aseq;
	$("#F1m2aseq").val(mhc2aseq);
	var mhc2bseq = document.getElementById("F1mhc2b").value;
	document.getElementById("F1m2bseq").value = mhc2bseq;
	$("#F1m2bseq").val(mhc2bseq);
    }


$("#loadbtnmhc1").click(function () {
    var mhc1aseq = document.getElementById("mhc1a").value;
    document.getElementById("m1aseq").value = mhc1aseq;
    $("#m1aseq").val(mhc1aseq);
});

$("#loadbtnmhc2").click(function () {
    var mhc2aseq = document.getElementById("mhc2a").value;
    document.getElementById("m2aseq").value = mhc2aseq;
    $("#m2aseq").val(mhc2aseq);
    var mhc2bseq = document.getElementById("mhc2b").value;
    document.getElementById("m2bseq").value = mhc2bseq;
    $("#m2bseq").val(mhc2bseq);
});

$("#loadbtn1").click(function () {
    var a_trav = document.getElementById("trav").value;
    var a_traj = document.getElementById("traj").value;
    var a_cdr3 = (document.getElementById("tacdrseq").value).toUpperCase();
    var aseq = a_trav+a_cdr3+a_traj;
    var b_trbv = document.getElementById("trbv").value;
    var b_trbj = document.getElementById("trbj").value;
    var b_cdr3 = (document.getElementById("tbcdrseq").value).toUpperCase();
    var bseq = b_trbv+b_cdr3+b_trbj;
    document.getElementById("taseq").value = aseq;
    document.getElementById("tbseq").value = bseq;
    $("#taseq").val(aseq);
    $("#tbseq").val(bseq);
    var mhc1aseq = document.getElementById("mhc1a").value;
    document.getElementById("m1aseq").value = mhc1aseq;
    $("#m1aseq").val(mhc1aseq);
    var pepseq = (document.getElementById("pseq").value).toUpperCase();
    document.getElementById("pseq").value = pepseq;
    $("#pseq").val(pepseq);
});

$("#loadbtn2").click(function () {
    var a_trav = document.getElementById("trav").value;
    var a_traj = document.getElementById("traj").value;
    var a_cdr3 = (document.getElementById("tacdrseq").value).toUpperCase();
    var aseq = a_trav+a_cdr3+a_traj;
    var b_trbv = document.getElementById("trbv").value;
    var b_trbj = document.getElementById("trbj").value;
    var b_cdr3 = (document.getElementById("tbcdrseq").value).toUpperCase();
    var bseq = b_trbv+b_cdr3+b_trbj;
    document.getElementById("taseq").value = aseq;
    document.getElementById("tbseq").value = bseq;
    $("#taseq").val(aseq);
    $("#tbseq").val(bseq);
    var mhc2aseq = document.getElementById("mhc2a").value;
    document.getElementById("m2aseq").value = mhc2aseq;
    $("#m2aseq").val(mhc2aseq);
    var mhc2bseq = document.getElementById("mhc2b").value;
    document.getElementById("m2bseq").value = mhc2bseq;
    $("#m2bseq").val(mhc2bseq);
    var pepseq = (document.getElementById("pseq").value).toUpperCase();
    document.getElementById("pseq").value = pepseq;
    $("#pseq").val(pepseq);
});

$("#loadbtntcra").click(function () {
    var a_trav = document.getElementById("trav").value;
    var a_traj = document.getElementById("traj").value;
    var a_cdr3 = (document.getElementById("tacdrseq").value).toUpperCase();
    var aseq = a_trav+a_cdr3+a_traj;
    document.getElementById("taseq").value = aseq;
    document.getElementById("trav_id").innerHTML = "";
    if( !a_cdr3 ) {
	document.getElementById("trav_id").innerHTML = "&#9888; CDR3 sequence not entered!";
    }
    if( !a_traj ) {
	document.getElementById("trav_id").innerHTML = "&#9888; TRAJ not selected!";
    }
    if( !a_trav ) {
	document.getElementById("trav_id").innerHTML = "&#9888; TRAV not selected!";
    }
});

$("#F2load_tcra").click(function () {
    var a_trav = document.getElementById("trav").value;
    var a_traj = document.getElementById("traj").value;
    var a_cdr3 = (document.getElementById("tacdrseq").value).toUpperCase();
    var aseq = a_trav+a_cdr3+a_traj;
    document.getElementById("taseq").value = aseq;
    document.getElementById("F2_trav_id").innerHTML = "";
    if( !a_cdr3 ) {
	document.getElementById("F2_trav_id").innerHTML = "&#9888; CDR3 sequence not entered!";
    }
    if( !a_traj ) {
	document.getElementById("F2_trav_id").innerHTML = "&#9888; TRAJ not selected!";
    }
    if( !a_trav ) {
	document.getElementById("F2_trav_id").innerHTML = "&#9888; TRAV not selected!";
    }
});

$("#F2load_tcrb").click(function () {
    var b_trbv = document.getElementById("trbv").value;
    var b_trbj = document.getElementById("trbj").value;
    var b_cdr3 = (document.getElementById("tbcdrseq").value).toUpperCase();
    var bseq = b_trbv+b_cdr3+b_trbj;
    document.getElementById("tbseq").value = bseq;
    document.getElementById("F2_trbv_id").innerHTML = "";
    if( !b_cdr3 ) {
	document.getElementById("F2_trbv_id").innerHTML = "&#9888; CDR3 sequence not entered!";
    }
    if( !b_trbj ) {
	document.getElementById("F2_trbv_id").innerHTML = "&#9888; TRBJ not selected!";
    }
    if( !b_trbv ) {
	document.getElementById("F2_trbv_id").innerHTML = "&#9888; TRBV not selected!";
    }
});

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

$("#loadbtntcrb").click(function () {
    var b_trbv = document.getElementById("trbv").value;
    var b_trbj = document.getElementById("trbj").value;
    var b_cdr3 = (document.getElementById("tbcdrseq").value).toUpperCase();
    var bseq = b_trbv+b_cdr3+b_trbj;
    document.getElementById("tbseq").value = bseq;
    document.getElementById("trbv_id").innerHTML = "";
    if( !b_cdr3 ) {
	document.getElementById("trbv_id").innerHTML = "&#9888; CDR3 sequence not entered!";
    }
    if( !b_trbj ) {
	document.getElementById("trbv_id").innerHTML = "&#9888; TRBJ not selected!";
    }
    if( !b_trbv ) {
	document.getElementById("trbv_id").innerHTML = "&#9888; TRBV not selected!";
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
    
$("#form1btn1").click(function () {
    $("#tachain").val("KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSWGKLQFGAGTQVVVTP");
    $("#tbchain").val("NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTE");
    $("#pchain").val("LLFGYPVYV");
    $("#m1achain").val("GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQR");
    $("#tachain").keyup();
    $("#tbchain").keyup();
    $("#pchain").keyup();
    $("#m1achain").keyup();
    //clear other forms
    $("#m2achain").val("");
    $("#m2achain").keyup();
    $("#m2bchain").val("");
    $("#m2bchain").keyup();
});
    
$("#form1btn2").click(function () {
    $("#tachain").val("MDAKTTQPNSMESNEEEPVHLPCNHSTISGTDYIHWYRQLPSQGPEYVIHGLTSNVNNRMASLAIAEDRKSSTLILHRATLRDAAVYYCILRDGRGGADGLTFGKGTHLIIQPYIQNP");
    $("#tbchain").val("DSGVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQGLQFLIQYYNGEERAKGNILERFSAQQFPDLHSELNLSSLELGDSALYFCASSVAVSAGTYEQYFGPGTRLTVTEDLKNVFP");
    $("#pchain").val("QQYPSGEGSFQPSQENPQ");
    $("#m2achain").val("EDIVADHVASYGVNLYQSYGPSGQYSHEFDGDEEFYVDLERKETVWQLPLFRRFRRFDPQFALTNIAVLKHNLNIVIKRSNSTAATNEVPEVTVFSKS");
    $("#m2bchain").val("SPEDFVYQFKGMCYFTNGTERVRLVTRYIYNREEYARFDSDVGVYRAVTPLGPPAAEYWNSQKEVLERTRAELDTVCRHNYQLELRTTLQRRVEPTVT");
    $("#tachain").keyup();
    $("#tbchain").keyup();
    $("#pchain").keyup();
    $("#m2achain").keyup();
    $("#m2bchain").keyup();
    //clear other forms
    $("#m1achain").val("");
    $("#m1achain").keyup();
});

$("#F1ExC1").click(function () {
    $("#tachain").val("KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSWGKLQFGAGTQVVVTP");
    $("#tbchain").val("NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTE");
    $("#pchain").val("LLFGYPVYV");
    $('select#F1mhc1a option:contains("HLA-A*02:01")').prop('selected',true);
    F1MHC1GeneSelection();
    //clear other forms
    $("#F1m2aseq").val("");
    $("#F1m2aseq").keyup();
    $("#F1m2bseq").val("");
    $("#F1m2bseq").keyup();
});

$("#F2ExC1").click(function () {
    $("#aspecies").prop('selectedIndex', 1);
    $("#tacdrseq").val("CAVGGSQGNLIF");
    $('select#trav option:contains("TRAV12-2*02")').prop('selected',true);
    $('select#traj option:contains("TRAJ1*01")').prop('selected',true);
    $("#F2load_tcra").click();
    $("#bspecies").prop('selectedIndex', 1);
    $("#tbcdrseq").val("CASSIRSSYEQYF");
    $('select#trbv option:contains("TRBV6-5*01")').prop('selected',true);
    $('select#trbj option:contains("TRBJ1-1*01")').prop('selected',true);
    $("#F2load_tcrb").click();
    $("#pseq").val("LLFGYPVYV");
    $('select#mhc1a option:contains("HLA-A*02:01")').prop('selected',true);
    MHC1GeneSelection();
    //clear other forms
    $("#m2aseq").val("");
    $("#m2aseq").keyup();
    $("#m2bseq").val("");
    $("#m2bseq").keyup();
});

$("#F1ExC2").click(function () {
    $("#tachain").val("MDAKTTQPNSMESNEEEPVHLPCNHSTISGTDYIHWYRQLPSQGPEYVIHGLTSNVNNRMASLAIAEDRKSSTLILHRATLRDAAVYYCILRDGRGGADGLTFGKGTHLIIQPYIQNP");
    $("#tbchain").val("DSGVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQGLQFLIQYYNGEERAKGNILERFSAQQFPDLHSELNLSSLELGDSALYFCASSVAVSAGTYEQYFGPGTRLTVTEDLKNVFP");
    $("#pchain").val("QQYPSGEGSFQPSQENPQ");
    $('select#F1mhc2a option:contains("HLA-DQA1*03:01")').prop('selected',true);
    $('select#F1mhc2b option:contains("HLA-DQB1*03:02")').prop('selected',true);
    F1MHC2GeneSelection();
    //clear other forms
    $("#F1m1aseq").val("");
    $("#F1m1aseq").keyup();
});


$("#F2ExC2").click(function () {
    $("#aspecies").prop('selectedIndex', 1);
    $("#tacdrseq").val("CILRDGRGGADGLTF");
    $('select#trav option:contains("TRAV26-2*01")').prop('selected',true);
    $('select#traj option:contains("TRAJ45*01")').prop('selected',true);
    $("#F2load_tcra").click();
    $("#bspecies").prop('selectedIndex', 1);
    $("#tbcdrseq").val("CASSVAVSAGTYEQYF");
    $('select#trbv option:contains("TRBV9*01")').prop('selected',true);
    $('select#trbj option:contains("TRBJ2-7*01")').prop('selected',true);
    $("#F2load_tcrb").click();
    $("#pseq").val("QQYPSGEGSFQPSQENPQ");
    $('select#mhc2a option:contains("HLA-DQA1*03:01")').prop('selected',true);
    $('select#mhc2b option:contains("HLA-DQB1*03:02")').prop('selected',true);
    MHC2GeneSelection();
    //clear other forms
    $("#m1aseq").val("");
    $("#m1aseq").keyup();
});


$("#ExBtn1").click(function () {
    $("#aspecies").prop('selectedIndex', 1);
    $("#bspecies").prop('selectedIndex', 1);
    $("#tacdrseq").val("CAVGGSQGNLIF");
    $("#tbcdrseq").val("CASSIRSSYEQYF");
    $('select#trav option:contains("TRAV12-2*02")').prop('selected',true);
    $('select#traj option:contains("TRAJ1*01")').prop('selected',true);
    $('select#trbv option:contains("TRBV6-5*01")').prop('selected',true);
    $('select#trbj option:contains("TRBJ1-1*01")').prop('selected',true);
    $("#pseq").val("LLFGYPVYV");
    $('select#mhc1a option:contains("HLA-A*02:01")').prop('selected',true);
    $("#loadbtntcra").click();
    $("#loadbtntcrb").click();
    MHC1GeneSelection();
    //clear other forms
    $("#m2aseq").val("");
    $("#m2aseq").keyup();
    $("#m2bseq").val("");
    $("#m2bseq").keyup();
});

$("#ExBtn2").click(function () {
    $("#aspecies").prop('selectedIndex', 1);
    $("#bspecies").prop('selectedIndex', 1);
    $("#tacdrseq").val("CILRDGRGGADGLTF");
    $("#tbcdrseq").val("CASSVAVSAGTYEQYF");
    $('select#trav option:contains("TRAV26-2*01")').prop('selected',true);
    $('select#traj option:contains("TRAJ45*01")').prop('selected',true);
    $('select#trbv option:contains("TRBV9*01")').prop('selected',true);
    $('select#trbj option:contains("TRBJ2-7*01")').prop('selected',true);
    $("#pseq").val("QQYPSGEGSFQPSQENPQ");
    $('select#mhc2a option:contains("HLA-DQA1*03:01")').prop('selected',true);
    $('select#mhc2b option:contains("HLA-DQB1*03:02")').prop('selected',true);
    $("#loadbtntcra").click();
    $("#loadbtntcrb").click();
    MHC2GeneSelection();
    //clear other forms
    $("#m1aseq").val("");
    $("#m1aseq").keyup();
});

$("#btn1").click(function () {
    $("#achain").val("KTTQPISMDSYEGQEVNITCSHNNIATNDYITWYQQFPSQGPRFIIQGYKTKVTNEVASLFIPADRKSSTLSLPRVSLSDTAVYYCLVGDMDQAGTALIFGKGTTLSVS");
    $("#bchain").val("VTQSPTHLIKTRGQQVTLRCSPKSGHDTVSWYQQALGQGPQFIFQYYEEEERQRGNFPDRFSGHQFPNYSSELNVNALLLGDSALYLCASSLGQTNYGYTFGSGTRLTVV");
    $("#achain").keyup();
    $("#bchain").keyup();
});
    
$("#btn2").click(function () {
    $("#sele_aspecies").prop('selectedIndex', 1);
    $("#sele_bspecies").prop('selectedIndex', 1);
    $("#sele_tacdrseq").val("CAVGGSQGNLIF");
    $("#sele_tbcdrseq").val("CASSIRSSYEQYF");
    $('select#sele_trav option:contains("TRAV8-6*01")').prop('selected',true);
    $('select#sele_traj option:contains("TRAJ42*01")').prop('selected',true);
    $('select#sele_trbv option:contains("TRBV19*01")').prop('selected',true);
    $('select#sele_trbj option:contains("TRBJ2-7*01")').prop('selected',true);
    $("#sele_loadbtntcra").click();
    $("#sele_loadbtntcrb").click();
});
    


$("#form2btn1").click(function () {
    $("#aspecies").prop('selectedIndex', 1);
    $("#bspecies").prop('selectedIndex', 1);
    $("#tacdrseq").val("CAVGGSQGNLIF");
    $("#tbcdrseq").val("CASSIRSSYEQYF");
    $('select#trav option:contains("TRAV12-2*02")').prop('selected',true);
    $('select#traj option:contains("TRAJ1*01")').prop('selected',true);
    $('select#trbv option:contains("TRBV6-5*01")').prop('selected',true);
    $('select#trbj option:contains("TRBJ1*01")').prop('selected',true);
    $("#pseq").val("LLFGYPVYV");
    $('select#mhc1a option:contains("HLA-A*02:01")').prop('selected',true);
    $("#loadbtn1").click();
    //clear other forms
    $("#m2aseq").val("");
    $("#m2aseq").keyup();
    $("#m2bseq").val("");
    $("#m2bseq").keyup();
});


$("#form2btn2").click(function () {
    $("#aspecies").prop('selectedIndex', 1);
    $("#bspecies").prop('selectedIndex', 1);
    $("#tacdrseq").val("CILRDGRGGADGLTF");
    $("#tbcdrseq").val("CASSVAVSAGTYEQYF");
    $('select#trav option:contains("TRAV26-2*01")').prop('selected',true);
    $('select#traj option:contains("TRAJ45*01")').prop('selected',true);
    $('select#trbv option:contains("TRBV9*01")').prop('selected',true);
    $('select#trbj option:contains("TRBJ2-7*01")').prop('selected',true);
    $("#pseq").val("QQYPSGEGSFQPSQENPQ");
    $('select#mhc2a option:contains("HLA-DQA1*03:01")').prop('selected',true);
    $('select#mhc2b option:contains("HLA-DQB1*03:02")').prop('selected',true);
    $("#loadbtn2").click();
    //clear other forms
    $("#m1aseq").val("");
    $("#m1aseq").keyup();
});


$("#tachain").keyup(function() {
	var x = document.getElementById("tachain").value;
	x = x.replace(/(\r\n|\n|\s|\r)/gm,"");
	var found = 0;
        travfile = $SCRIPT_ROOT + '/static/genemapdata/TRAV_HUMAN.json';
        $.getJSON(travfile, function (data) {
            $.each(data, function (key, entry) {
                    var re = new RegExp(entry.seq.substring(3))
		    if (re.test(x))
			{
				document.getElementById("trav_id").innerHTML = "Human " + entry.name + " gene identified";
				found = 1;
				return true;
			}
                })
                });		
	if (found == 0) {  document.getElementById("trav_id").innerHTML = ""; }
	});
    
$("#tbchain").keyup(function() {
        var x = document.getElementById("tbchain").value;
        x = x.replace(/(\r\n|\n|\s|\r)/gm,"");
	var found = 0;
        travfile = $SCRIPT_ROOT + '/static/genemapdata/TRBV_HUMAN.json';
        $.getJSON(trbvfile, function (data) {
            $.each(data, function (key, entry) {
                    var re = new RegExp(entry.seq.substring(3))
                    if (re.test(x))
                        {
                                document.getElementById("trbv_id").innerHTML = "Human " + entry.name + " gene identified";
                                found = 1;
                                return true;
                        }
                })
                });
        if (found == 0) {  document.getElementById("trbv_id").innerHTML = ""; }
        });

$("#reset1").click(function () {
	document.getElementById("trav_id").innerHTML = "";
	document.getElementById("trbv_id").innerHTML = "";
	});




});

