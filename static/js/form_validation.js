$(document).ready(function() {

    $("#submit1").click(function() {
        isValidSeq("tachain");
        isValidSeq("tbchain");
        checkPep("pchain", "I", 8);
        checkMhc1("F1mhc1speciestype", "F1mhc1a", "F1mhc2a", "F1mhc2b", "F1m1aseq");
    });

    $("#submit2").click(function() {
        isValidSeq("tachain");
        isValidSeq("tbchain");
        checkPep("pchain", "II", 9);
        checkMhc2("F1mhc2speciestype", "F1mhc1a", "F1mhc2a", "F1mhc2b", "F1m2aseq", "F1m2bseq");
    });

    $("#submit3").click(function() {
        isValidSeq("tacdrseq");
        isValidSeq("tbcdrseq");
        checkPep("pseq", "I", 8);
        checkMhc1("mhc1speciestype", "mhc1a", "mhc2a", "mhc2b", "m1aseq");
    });

    $("#submit4").click(function() {
        isValidSeq("tacdrseq");
        isValidSeq("tbcdrseq");
        checkPep("pseq", "II", 9);
        checkMhc2("mhc2speciestype", "mhc1a", "mhc2a", "mhc2b", "m2aseq", "m2bseq");
    });

    $("#submit5").click(function() {
        isValidSeq("achain");
        isValidSeq("bchain");
    });

    $("#submit6").click(function() {
        isValidSeq("sele_tacdrseq");
        isValidSeq("sele_tbcdrseq");
    });

    function isValidSeq(id) {
        var field = document.getElementById(id);
        if (!field.value) {
            field.setCustomValidity("Please enter the sequence.");
            return false;
        }
        if (!/^[ARNDCEQGHILKMFPSTWYV]+$/.exec(field.value)) {
            field.setCustomValidity("Sequence can only contain standard amino acids.");
            return false;
        }
        field.setCustomValidity("");
        return true;
    }

    function checkPep(id, mhcClass, min) {
        if (isValidSeq(id)) {
            var pep = document.getElementById(id);
            if (pep.value.length < min) {
                var msg = "Currently MHC " + mhcClass + " complex modeling supports only " + min + "-mer peptides.";
                pep.setCustomValidity(msg);
            }
        }
    }

    function isSpeciesPicked(id) {
        var species = document.getElementById(id);
        if (!species.value) {
            species.setCustomValidity("Please select the species.");
            return false;
        }
        species.setCustomValidity("");
        return true;
    } 

    function setMhcErrorMsgs(fieldIds, shouldSetMsgs) {
        for (var i = 0; i < fieldIds.length; i++) {
            if (shouldSetMsgs[i])
                document.getElementById(fieldIds[i]).setCustomValidity("Please select or enter the MHC sequence(s).");
            else
                document.getElementById(fieldIds[i]).setCustomValidity("");
        }
    }

    function checkMhc1(species, menu1, menu2a, menu2b, id) {
        if (!isSpeciesPicked(species))
            return;
        var mhc1 = document.getElementById(id);
        if (!mhc1.value)
            setMhcErrorMsgs([menu1, menu2a, menu2b], [true, false, false]);
        else
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, false, false]);
    }

    function checkMhc2(species, menu1, menu2a, menu2b, id2a, id2b) {
        if (!isSpeciesPicked(species))
            return;
        var mhc2a = document.getElementById(id2a);
        var mhc2b = document.getElementById(id2b);
        if (!mhc2a.value && !mhc2b.value)
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, true, true]);
        else if (!mhc2a.value)
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, true, false]);
        else if (!mhc2b.value)
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, false, true]);
        else
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, false, false]);
    }
});
