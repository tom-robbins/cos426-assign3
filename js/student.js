var Student = {
    // please fill in your name and NetID
    // your NetID is the part of your email before @princeton.edu
    'name'  : 'Porter Sherman and Tom Robbins',
    'netID' : 'wps3 thomasrr',
};

Student.updateHTML = function( ) {
    var studentInfo = this.name + ' &lt;' + this.netID + '&gt;';
    document.getElementById('student').innerHTML = studentInfo;
}
