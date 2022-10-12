var parent = document.getElementById('vwr');
var viewer = pv.Viewer(parent, 
		       { quality : 'medium', width: 'auto', height : 'auto',
			 antialias : true, outline : true});
var structure;
var symm_axis;

function lines() {
    viewer.clear();
    var tcr = structure.select({cnames: 'AB'});
    viewer.cartoon('tcr', tcr, { color: pv.color.byChain() } );
    viewer.lines('tcr', tcr, { color: pv.color.byChain() });
}

function view_lines() {
    viewer.clear();
    var tcr = structure.select({cnames: 'AB'});
    viewer.cartoon('tcr', tcr, { color: pv.color.byChain() } );
    viewer.lines('tcr', tcr, { color: pv.color.byChain() });
}

function cdr3lines() {
    viewer.clear();	
    var tcr = structure.select({cnames: 'AB'});
    viewer.cartoon('tcr', tcr, { color: pv.color.byChain() });
    var cdr3 = structure.select({ rnumRange : ['107','139'] });
    viewer.lines('cdr3', cdr3, { color: pv.color.byElement() });
}

function view_cdr3lines() {
    viewer.clear();	
    var tcr = structure.select({cnames: 'AB'});
    viewer.cartoon('tcr', tcr, { color: pv.color.byChain() });
    var cdr3 = structure.select({ rnumRange : ['107','139'] });
    viewer.lines('cdr3', cdr3, { color: pv.color.byElement() });
}

function cartoon() {
    viewer.clear();
    var tcr = structure.select({cnames: 'ABCDE'});
    viewer.cartoon('tcr', tcr, { color: pv.color.byChain() } );
}

function reset() { 
    viewer.clear();
    cartoon();
    viewer.autoZoom();
}

function tcr() {
    var xhr = new XMLHttpRequest();
    xhr.open('GET', modelfname );
    xhr.onreadystatechange = function() {
	if (xhr.readyState == 4) {
	    structure = pv.io.pdb(xhr.responseText);
	    cartoon();
	    viewer.autoZoom();
	}
    }
    xhr.send();
    document.getElementById('pickedatom').innerHTML = 'TCR model loaded';
}

function load_tmplt(tmpltfname) {
    pv.io.fetchPdb(tmpltfname, function(structure) {
            selep = structure.select('protein');
            viewer.on('viewerReady', function() {
                    viewer.cartoon('selep',selep , { color: pv.color.rainbow() } );
                });

        });
}


function read_tmplt_Acdr1() {
    var tmpltfname = tj.rundir_spath+"/"+jobid+"/"+tj.acdr1_tmplt_pdb+"_Acdr1_tmplt.pdb";
    pv.io.fetchPdb(tmpltfname, function(structure) {
            var acdr1 = structure.select({ cnames: tj.acdr1_tmplt_pdb_chain, rnumRange : ['24','43'] });
	    viewer.on('viewerReady', function() {
		    viewer.cartoon('acdr1', acdr1, { color: pv.color.uniform('green')});
		});
	    
	});
}
function read_tmplt_Acdr2() {
    var tmpltfname = tj.rundir_spath+"/"+jobid+"/"+tj.acdr2hv4_tmplt_pdb+"_Acdr2hv4_tmplt.pdb";
    pv.io.fetchPdb(tmpltfname, function(structure) {
            var acdr2 = structure.select({ cnames: tj.acdr2hv4_tmplt_pdb_chain, rnumRange : ['56','91'] });
	    viewer.on('viewerReady', function() {
		    viewer.cartoon('acdr2', acdr2, { color: pv.color.uniform('green')});
		});
	    
	});
}
function read_tmplt_Acdr3() {
    var tmpltfname = tj.rundir_spath+"/"+jobid+"/"+tj.acdr3_tmplt_pdb+"_Acdr3_tmplt.pdb";
    pv.io.fetchPdb(tmpltfname, function(structure) {
            var acdr3 = structure.select({ cnames: tj.acdr3_tmplt_pdb_chain, rnumRange : ['107','139'] });
	    viewer.on('viewerReady', function() {
		    viewer.cartoon('acdr3', acdr3, { color: pv.color.uniform('green')});
		});
	    
	});
}
function read_tmplt_Bcdr1() {
    var tmpltfname = tj.rundir_spath+"/"+jobid+"/"+tj.bcdr1_tmplt_pdb+"_Bcdr1_tmplt.pdb";
    pv.io.fetchPdb(tmpltfname, function(structure) {
            var bcdr1 = structure.select({ cnames: tj.bcdr1_tmplt_pdb_chain, rnumRange : ['24','43'] });
	    viewer.on('viewerReady', function() {
		    viewer.cartoon('bcdr1', bcdr1, { color: pv.color.uniform('green')});
		});
	    
	});
}
function read_tmplt_Bcdr2() {
    var tmpltfname = tj.rundir_spath+"/"+jobid+"/"+tj.bcdr2hv4_tmplt_pdb+"_Bcdr2hv4_tmplt.pdb";
    pv.io.fetchPdb(tmpltfname, function(structure) {
            var bcdr2 = structure.select({ cnames: tj.bcdr2hv4_tmplt_pdb_chain, rnumRange : ['56','91'] });
	    viewer.on('viewerReady', function() {
		    viewer.cartoon('bcdr2', bcdr2, { color: pv.color.uniform('green')});
		});
	    
	});
}
function read_tmplt_Bcdr3() {
    var tmpltfname = tj.rundir_spath+"/"+jobid+"/"+tj.bcdr3_tmplt_pdb+"_Bcdr3_tmplt.pdb";
    pv.io.fetchPdb(tmpltfname, function(structure) {
            var bcdr3 = structure.select({ cnames: tj.bcdr3_tmplt_pdb_chain, rnumRange : ['107','139'] });
	    viewer.on('viewerReady', function() {
		    viewer.cartoon('bcdr3', bcdr3, { color: pv.color.uniform('green')});
		});
	    
	});
}
function read_tmplt_ori() {
    var tmpltfnameA = tj.rundir_spath+"/"+jobid+"/"+tj.ori_tmplt_Apdb+"_oriA_tmplt.pdb";
    load_tmplt(tmpltfnameA);
    var tmpltfnameB = tj.rundir_spath+"/"+jobid+"/"+tj.ori_tmplt_Bpdb+"_oriB_tmplt.pdb";
    load_tmplt(tmpltfnameB);
}



document.getElementById('cartoon').onclick = cartoon;
//document.getElementById('lines').onclick = lines;
//document.getElementById('cdr3lines').onclick = cdr3lines;  

if (document.getElementById('reset') != null) {
    document.getElementById('reset').onclick = reset;
}


//if (document.getElementById('tmplt_Acdr1') != null) {
    //document.getElementById('tmplt_Acdr1').onclick = read_tmplt_Acdr1;
//}

//if (document.getElementById('tmplt_Acdr2') != null) {
    //document.getElementById('tmplt_Acdr2').onclick = read_tmplt_Acdr2;
//}
//if (document.getElementById('tmplt_Acdr3') != null) {
    //document.getElementById('tmplt_Acdr3').onclick = read_tmplt_Acdr3;
//}
//if (document.getElementById('tmplt_Bcdr1') != null) {
    //document.getElementById('tmplt_Bcdr1').onclick = read_tmplt_Bcdr1;
//}
//if (document.getElementById('tmplt_Bcdr2') != null) {
//  document.getElementById('tmplt_Bcdr2').onclick = read_tmplt_Bcdr2;
//}
//if (document.getElementById('tmplt_Bcdr3') != null) {
//  document.getElementById('tmplt_Bcdr3').onclick = read_tmplt_Bcdr3;
//}
//if (document.getElementById('tmplt_ori') != null) {
//  document.getElementById('tmplt_ori').onclick = read_tmplt_ori;
//}



function setColorForAtom(go, atom, color) {
    var view = go.structure().createEmptyView();
    view.addAtom(atom);
    go.colorBy(pv.color.uniform(color), view);
}

// variable to store the previously picked atom. Required for resetting the color
// whenever the mouse moves.
var prevPicked = null;
// add mouse move event listener to the div element containing the viewer. Whenever
// the mouse moves, use viewer.pick() to get the current atom under the cursor.
parent.addEventListener('mousemove', function(event) {
	var rect = viewer.boundingClientRect();
	var picked = viewer.pick({ x : event.clientX - rect.left,
				   y : event.clientY - rect.top });
	if (prevPicked !== null && picked !== null &&
	    picked.target() === prevPicked.atom) {
	    return;
	}
	if ((picked.target().qualifiedName().substring(0, 1) === "G") ||
	    (picked.target().qualifiedName().substring(0, 1) === "H")) { return; }
	if (prevPicked !== null) {
	    // reset color of previously picked atom.
	    setColorForAtom(prevPicked.node, prevPicked.atom, prevPicked.color);
	}
	if (picked !== null) {
	    var atom = picked.target();
	    document.getElementById('pickedatom').innerHTML = atom.qualifiedName();
	    // get RGBA color and store in the color array, so we know what it was
	    // before changing it to the highlight color.
	    var color = [0,0,0,0];
	    picked.node().getColorForAtom(atom, color);
	    prevPicked = { atom : atom, color : color, node : picked.node() };
	    
	    setColorForAtom(picked.node(), atom, 'magenta');
	} else {
	    document.getElementById('pickedatom').innerHTML = '&nbsp;';
	    prevPicked = null;
	}
	viewer.requestRedraw();
    });

window.onresize = function(event) {
    viewer.fitParent();
}
    document.addEventListener('DOMContentLoaded', tcr);


$(document).ready(function() {

	function show_afw(){
	    viewer.clear();
	    var afw1 = structure.select({ cname : 'A', rnumRange : ['1','24'] });
	    var afw2 = structure.select({ cname : 'A', rnumRange : ['43','56'] });
	    var afw3 = structure.select({ cname : 'A', rnumRange : ['91','107'] });
	    var afw4 = structure.select({ cname : 'A', rnumRange : ['139','149'] });
	    viewer.cartoon('afw1', afw1, { color: pv.color.uniform('red')});
	    viewer.cartoon('afw2', afw2, { color: pv.color.uniform('red')});
	    viewer.cartoon('afw3', afw3, { color: pv.color.uniform('red')});
	    viewer.cartoon('afw4', afw4, { color: pv.color.uniform('red')});
	}
	function show_acdr1(){		
	    viewer.clear();
	    var acdr1 = structure.select({ cname : 'A', rnumRange : ['24','43'] });
	    viewer.cartoon('acdr1', acdr1, { color: pv.color.uniform('red')});
	}
	function show_acdr2hv4(){		
	    viewer.clear();
	    var acdr2hv4 = structure.select({ cname : 'A', rnumRange : ['56','91'] });
	    viewer.cartoon('acdr2hv4', acdr2hv4, { color: pv.color.uniform('red')});
	}
	function show_acdr3(){		
	    viewer.clear();
	    var acdr1 = structure.select({ cname : 'A', rnumRange : ['107','139'] });
	    viewer.cartoon('acdr1', acdr1, { color: pv.color.uniform('red')});
	}
	function show_bfw(){
	    viewer.clear();
	    var bfw1 = structure.select({ cname : 'B', rnumRange : ['1','24'] });
	    var bfw2 = structure.select({ cname : 'B', rnumRange : ['43','56'] });
	    var bfw3 = structure.select({ cname : 'B', rnumRange : ['91','107'] });
	    var bfw4 = structure.select({ cname : 'B', rnumRange : ['139','149'] });
	    viewer.cartoon('bfw', bfw1, { color: pv.color.uniform('blue')});
	    viewer.cartoon('bfw', bfw2, { color: pv.color.uniform('blue')});
	    viewer.cartoon('bfw', bfw3, { color: pv.color.uniform('blue')});
	    viewer.cartoon('bfw', bfw4, { color: pv.color.uniform('blue')});
	}

	function show_bcdr1(){		
	    viewer.clear();
	    var bcdr1 = structure.select({ cname : 'B', rnumRange : ['24','43'] });
	    viewer.cartoon('bcdr1', bcdr1, { color: pv.color.uniform('blue')});
	}
	function show_bcdr2hv4(){		
	    viewer.clear();
	    var bcdr2hv4 = structure.select({ cname : 'B', rnumRange : ['56','91'] });
	    viewer.cartoon('bcdr2hv4', bcdr2hv4, { color: pv.color.uniform('blue')});
	}
	function show_bcdr3(){		
	    viewer.clear();
	    var bcdr3 = structure.select({ cname : 'B', rnumRange : ['107','139'] });
	    viewer.cartoon('bcdr3', bcdr3, { color: pv.color.uniform('blue')});
	}

	var isClicked = false;
	$('#afw, #acdr1, #acdr2hv4, #acdr3, #bfw, #bcdr1, #bcdr2hv4, #bcdr3').click(function() {
		if ( ($(this).hasClass("clicked-once")) && (isClicked) ) {
		    cartoon();
		    $('#afw, #acdr1, #acdr2hv4, #acdr3, #bfw, #bcdr1, #bcdr2hv4, #bcdr3').removeClass("clicked-once");
		    isClicked = false;
		    $('#afw, #acdr1, #acdr2hv4, #acdr3, #bfw, #bcdr1, #bcdr2hv4, #bcdr3').css('font-weight', 'normal');
		}
		else {
		    //alert($(this).attr("id"));
		    $('#afw, #acdr1, #acdr2hv4, #acdr3, #bfw, #bcdr1, #bcdr2hv4, #bcdr3').removeClass("clicked-once");
		    $(this).addClass("clicked-once");
		    isClicked = true;
		    var tcrseg = $(this).attr("id");
		    eval("show_" + tcrseg + "()");
		    $('#afw, #acdr1, #acdr2hv4, #acdr3, #bfw, #bcdr1, #bcdr2hv4, #bcdr3').css('font-weight', 'normal');
		    $(this).css('font-weight', 'bold');
		}
	    });
	
	$('#afw, #acdr1, #acdr2hv4, #acdr3, #bfw, #bcdr1, #bcdr2hv4, #bcdr3').on('mouseover',function() {
		if(!isClicked){
		    var tcrseg = $(this).attr("id");
		    eval("show_" + tcrseg + "()");
		    $(this).css('font-weight', 'bold');
		};
            });

	$('#afw, #acdr1, #acdr2hv4, #acdr3, #bfw, #bcdr1, #bcdr2hv4, #bcdr3').mouseout(function(){
		if(!isClicked) {
		    cartoon();
		    $(this).css('font-weight', 'normal');
		}
	    });

	$('#cdr3lines, #lines').click(function() {
		if ( ($(this).hasClass("clicked-once")) && (isClicked) ) {
		    reset();
		    $('#cdr3lines, #lines').removeClass("clicked-once");
		    isClicked = false;
		}
		else {
		    //alert($(this).attr("id"));
		    $('#cdr3lines, #lines').removeClass("clicked-once");
		    $(this).addClass("clicked-once");
		    isClicked = true;
		    var view_func = $(this).attr("id");
		    eval("view_" + view_func + "()");
		}
	    });


	$('#tmplt_Acdr1, #tmplt_Acdr2, #tmplt_Acdr3, #tmplt_Bcdr1, #tmplt_Bcdr2, #tmplt_Bcdr3, #tmplt_ori').click(function() {
		if ( ($(this).hasClass("clicked-once")) && (isClicked) ) {
		    //hide template function goes here
		    $('#tmplt_Acdr1, #tmplt_Acdr2, #tmplt_Acdr3, #tmplt_Bcdr1, #tmplt_Bcdr2, #tmplt_Bcdr3, #tmplt_ori').removeClass("clicked-once");
		    isClicked = false;
		}
		else {
		    //alert($(this).attr("id"));
		    $('#tmplt_Acdr1, #tmplt_Acdr2, #tmplt_Acdr3, #tmplt_Bcdr1, #tmplt_Bcdr2, #tmplt_Bcdr3, #tmplt_ori').removeClass("clicked-once");
		    $(this).addClass("clicked-once");
		    isClicked = true;
		    var view_func = $(this).attr("id");
		    eval("read_" + view_func + "()");
		}
	    });

    });


