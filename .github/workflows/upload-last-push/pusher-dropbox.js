const dropboxV2Api = require('dropbox-v2-api');
const path = require('path');
const fs = require('fs');
// git head hash is stored in the log
var git_hash = require('child_process').execSync('git rev-parse HEAD').toString();
// Token used to push files to dropbox
const token = process.env.DROPBOX_PUSH_KEY;
const dropbox = dropboxV2Api.authenticate({
    token: token
});

function UploadFile(filename, isLog = false) {

    var file_raw = path.basename(filename);
    if (isLog) file_raw = "logs/"+file_raw;
    const dropboxUploadStream = dropbox({
        resource: 'files/upload',
        parameters: {
            path: '/' + file_raw,
            mode: "overwrite"
        }
    }, (err, result, response) => {
        if (!isLog) {
            let log_filename = logs_dir + '/' + GetTimestamp() + ".log";
            let log_content = {};
            if (err === undefined || err === null) {
                // If success -> 
                result.git_hash = git_hash;
                log_content = result;
            } else {
                err.git_hash = git_hash;
                log_content = err;
            }
            console.log(log_content);
            let data = JSON.stringify(log_content);
            fs.writeFileSync(log_filename, data);
            UploadFile(log_filename, true);
        } else {
            if (err !== undefined && err !== null) {
                // If log storage fails :(
                console.log(err)
                return -1;
            }
        }
    });

    fs.createReadStream(filename).pipe(dropboxUploadStream);
}

// push to dropbox the delivery file
var compressed_kratos = '/zip/latest-linux-x64.tgz';
// Create logs folder
var logs_dir = '/zip/logs';

var isWin = process.platform === "win32";
if (isWin) {
    compressed_kratos = 'D:/a/latest-windows-x64.zip';
    logs_dir = 'D:/a/logs';
}

if (!fs.existsSync(logs_dir)) {
    fs.mkdirSync(logs_dir);
}
UploadFile(compressed_kratos);


function GetTimestamp() {
    var ts_hms = new Date();

    var total = ts_hms.getFullYear() + '-' + 
    ("0" + (ts_hms.getMonth() + 1)).slice(-2) + '-' + 
    ("0" + (ts_hms.getDate())).slice(-2) + '-' +
    ("0" + ts_hms.getHours()).slice(-2) + '-' +
    ("0" + ts_hms.getMinutes()).slice(-2) + '-' +
    ("0" + ts_hms.getSeconds()).slice(-2);

    return total;
}
return 0;




// # - name: Zip generate
// #   uses: montudor/action-zip@v0.1.0
// #   with:
// #     args: zip -qq -r ./latest.zip /zip/dist/runkratos