function openCity(evt, cityName) {
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
	tabcontent[i].style.display = "none";
    }
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
	tablinks[i].className = tablinks[i].className.replace(" active", "");
    }
    document.getElementById(cityName).style.display = "block";
    evt.currentTarget.className += " active";
}


$(document).ready(function(){


var travfile;
var trajfile;
var trbvfile;
var trbjfile;
AGeneSelection();
BGeneSelection();

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

document.getElementById("bspecies").onchange = function() {BGeneSelection()};
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
    trbvdropdown.append('<option selected="true" disabled value="">Choose TRBV</option>');
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



    $("#loadbtn").click(function () {
	var a_trav = document.getElementById("trav").value;
	var a_traj = document.getElementById("traj").value;
	var a_cdr3 = (document.getElementById("acdr").value).toUpperCase();
	var aseq = a_trav+a_cdr3+a_traj;
	var b_trbv = document.getElementById("trbv").value;
	var b_trbj = document.getElementById("trbj").value;
	var b_cdr3 = (document.getElementById("bcdr").value).toUpperCase();
	var bseq = b_trbv+b_cdr3+b_trbj;
	document.getElementById("alphachain").value = aseq;
	document.getElementById("betachain").value = bseq;
	$("#alphachain").val(aseq);
	$("#betachain").val(bseq);
    });

document.getElementById("defaultOpen").click();
$("#flip").click(function(){
	$("#panel").slideToggle("slow");
    });
$("#flip2").click(function(){
	$("#panel2").slideToggle("slow");
    });
$("#flip3").click(function(){
	$("#panel3").slideToggle("slow");
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
	$("#acdr").val("CAVGGSQGNLIF");
	$("#bcdr").val("CASSIRSSYEQYF");
        $('select#trav option:contains("TRAV8-6*01")').prop('selected',true);
        $('select#traj option:contains("TRAJ42*01")').prop('selected',true);
        $('select#trbv option:contains("TRBV19*01")').prop('selected',true);
        $('select#trbj option:contains("TRBJ2-7*01")').prop('selected',true);
	$("#loadbtn").click();
    });


$("#achain").keyup(function() {
    var x = document.getElementById("achain").value;
    x = x.replace(/(\r\n|\n|\s|\r)/gm,"");
    var found = 0;
    travfile = $SCRIPT_ROOT + '/static/genemapdata/TRAV_HUMAN.json';
    $.getJSON(travfile, function (data) {
        $.each(data, function (key, entry) {
            var re = new RegExp(entry.seq.substring(3))
	    if (re.test(x))
	    {
		document.getElementById("trav_id").innerHTML = "Human " + entry.name + " gene identified";
		document.getElementsByClassName("trav_id").innerHTML = "Human " + entry.name + " gene identified";
		found = 1;
		return true;
	    }
        })
            });		
    if (found == 0) {  
	document.getElementById("trav_id").innerHTML = "";
	document.getElementByClassName("trav_id").innerHTML = "";
    }
});
    
$("#bchain").keyup(function() {
        var x = document.getElementById("bchain").value;
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
	document.getElementByClassName("trav_id").innerHTML = "";
	document.getElementById("trbv_id").innerHTML = "";
	});

});

