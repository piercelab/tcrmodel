$(document).ready(function(){
    $("#loadbtn").click(function () {
        var a_trav = document.getElementById("trav").value;
        var a_traj = document.getElementById("traj").value;
        var a_cdr3 = (document.getElementById("acdr").value).toUpperCase();
        var aseq = a_trav+a_cdr3+a_traj;
        var b_trbv = document.getElementById("trbv").value;
        var b_trbj = document.getElementById("trbj").value;
        var b_cdr3 = (document.getElementById("bcdr").value).toUpperCase();
        var bseq = b_trbv+b_cdr3+b_trbj;
	alert(a_trav);
        document.getElementById("alphachain").value = aseq;
        document.getElementById("betachain").value = bseq;
        $("#alphachain").val(aseq);
        $("#betachain").val(bseq);
    });
    
    document.getElementById("opensecondtab").click();
  //  AGeneSelection();
  //  BGeneSelection();
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
