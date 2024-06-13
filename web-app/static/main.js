document.addEventListener('DOMContentLoaded', (event) => {
    const inputField = document.getElementById("dot_bracket");
    const errorMessage = document.getElementById("error_message");

    inputField.addEventListener('input', () => {
        errorMessage.textContent = ''; // Clear the error message
        inputField.style.borderColor = ''; // Reset border color

        // Filter out invalid characters
        const validCharacters = ['.', '(', ')', '[', ']', '{', '}'];
        let filteredValue = '';
        for (let char of inputField.value) {
            if (validCharacters.includes(char)) {
                filteredValue += char;
            }
        }
        inputField.value = filteredValue;
    });
});

function validateDotBracket() {
    const inputField = document.getElementById("dot_bracket");
    const errorMessage = document.getElementById("error_message");
    const input = inputField.value;
    const validCharacters = new Set(['.', '(', ')', '[', ']', '{', '}']);
    const stack = [];
    const pairs = {')': '(', ']': '[', '}': '{'};

    for (let char of input) {
        if (!validCharacters.has(char)) {
            errorMessage.textContent = `Invalid character found: ${char}`;
            inputField.style.borderColor = 'red';
            return false;
        }
        if (char === '(' || char === '[' || char === '{') {
            stack.push(char);
        } else if (char === ')' || char === ']' || char === '}') {
            if (stack.length === 0 || stack.pop() !== pairs[char]) {
                errorMessage.textContent = 'Unbalanced brackets detected!';
                inputField.style.borderColor = 'red';
                return false;
            }
        }
    }

    if (stack.length !== 0) {
        errorMessage.textContent = 'Unbalanced brackets detected!';
        inputField.style.borderColor = 'red';
        return false;
    }

    return true;
}

function downloadResult() {
    const resultDiv = document.querySelector('.result ol');
    if (!resultDiv) return;

    const structure = document.querySelector('.result h3:last-of-type').textContent.trim();
    const verifyDesign = document.querySelector('#verify_design').checked;
    let text = `${structure}\n`;
    const items = resultDiv.querySelectorAll('li');

    items.forEach((item, index) => {
        const sequenceText = item.querySelector('.design-sequence').textContent.trim();
        if (verifyDesign) {
            const isInvalid = sequenceText.includes('invalid');
            if (!isInvalid) {
                text += `${sequenceText.split(' ')[0].trim()}\n`;
            }
        } else {
            text += `${sequenceText}\n`;
        }
    });

    const blob = new Blob([text], { type: 'text/plain' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = 'rna_design_results.txt';
    link.click();
}


function clearResult() {
    window.location.href = '/';
}