
document.addEventListener('DOMContentLoaded', () => {
    document.querySelectorAll('.toggle-btn').forEach(btn => {
        btn.addEventListener('click', function() {
            const card = this.closest('.result-card');
            card.classList.toggle('expanded');

            // Change button text
            this.textContent = card.classList.contains('expanded') 
                ? 'Show less' 
                : 'Show more';
        });
    });
});



