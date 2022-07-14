$(document).ready(function() {

    $("#submit1").click(function() {
        isNotEmptyValid("tachain");
        isNotEmptyValid("tbchain");
        checkPep("pchain", "I", 8);
        checkMhc1("F1mhc1a", "F1mhc2a", "F1mhc2b", "F1m1aseq", "F1mhc1area");
    });

    $("#submit2").click(function() {
        isNotEmptyValid("tachain");
        isNotEmptyValid("tbchain");
        checkPep("pchain", "II", 9);
        checkMhc2("F1mhc1a", "F1mhc2a", "F1mhc2b", "F1m2aseq", "F1m2bseq", "F1mhc2area");
    });

    $("#submit3").click(function() {
        isNotEmptyValid("tacdrseq");
        isNotEmptyValid("tbcdrseq");
        checkPep("pseq", "I", 8);
        checkMhc1("mhc1a", "mhc2a", "mhc2b", "m1aseq", "mhc1area");
    });

    $("#submit4").click(function() {
        isNotEmptyValid("tacdrseq");
        isNotEmptyValid("tbcdrseq");
        checkPep("pseq", "II", 9);
        checkMhc2("mhc1a", "mhc2a", "mhc2b", "m2aseq", "m2bseq", "mhc2area");
    });

    $("#submit5").click(function() {
        isNotEmptyValid("achain");
        isNotEmptyValid("bchain");
    });

    $("#submit6").click(function() {
        isNotEmptyValid("sele_tacdrseq");
        isNotEmptyValid("sele_tbcdrseq");
    });

    function isValidSeq(field) {
        if (!/^[ARNDCEQGHILKMFPSTWYV]+$/.exec(field.value)) {
            field.setCustomValidity("Sequence can only contain standard amino acids.");
            return false;
        }
        field.setCustomValidity("");
        return true;
    }

    function isNotEmptyValid(id) {
        var field = document.getElementById(id);
        if (!field.value) {
            field.setCustomValidity("Please enter the sequence.");
            return false;
        }
        return isValidSeq(field);
    }

    function checkPep(id, mhcClass, min) {
        if (isNotEmptyValid(id)) {
            var pep = document.getElementById(id);
            if (pep.value.length < min) {
                var msg = "Currently MHC " + mhcClass + " complex modeling supports only " + min + "-mer peptides.";
                pep.setCustomValidity(msg);
            }
        }
    }

    function setMhcErrorMsgs(fieldIds, shouldSetMsgs) {
        for (var i = 0; i < fieldIds.length; i++) {
            if (shouldSetMsgs[i])
                document.getElementById(fieldIds[i]).setCustomValidity("Please select or enter the MHC sequence(s).");
            else
                document.getElementById(fieldIds[i]).setCustomValidity("");
        }
    }

    function checkMhc1(menu1, menu2a, menu2b, id, collapseArea) {
        var mhc1 = document.getElementById(id);
        if (!mhc1.value)
            setMhcErrorMsgs([menu1, menu2a, menu2b], [true, false, false]);
        else {
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, false, false]);
            if (!isValidSeq(mhc1))
                $("#"+collapseArea).prop("class", "px-0 collapse show");
        }
    }

    function checkMhc2(menu1, menu2a, menu2b, id2a, id2b, collapseArea) {
        var mhc2a = document.getElementById(id2a);
        var mhc2b = document.getElementById(id2b);
        if (!mhc2a.value && !mhc2b.value)
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, true, true]);
        else if (!mhc2a.value)
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, true, false]);
        else if (!mhc2b.value)
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, false, true]);
        else {
            setMhcErrorMsgs([menu1, menu2a, menu2b], [false, false, false]);
            if (!isValidSeq(mhc2a) || !isValidSeq(mhc2b))
                $("#"+collapseArea).prop("class", "px-0 collapse show");
        }
    }
});
